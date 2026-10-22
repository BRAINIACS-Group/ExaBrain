// node_nd_ordering.h
// Nested-Dissection-Ordering (METIS_NodeND) auf dem KNOTENGRAPHEN statt auf dem DoF-Graphen.
//
// Idee: Alle DoFs eines Gitterpunkts (u_x,u_y,u_z und ggf. p an Hex-Ecken) werden zu EINEM
// Vertex zusammengefasst (Gewicht = Anzahl DoFs). Der Graph hat bei Q2/Q1 etwa 1/3 der Vertices
// und 1/9 der Kanten (40 Mio DoFs: ca. 12.8 Mio Vertices, ca. 0.9e9 Eintraege statt 8.5e9), passt
// also in 32-bit-METIS (idx_t = int) und braucht auf Rank 0 nur ca. 20-40 GB statt 300-500 GB.
// Die Knotenpermutation wird wieder auf die DoFs expandiert (DoFs eines Knotens stehen
// hintereinander in der Pivotfolge) und als PERM_IN an MUMPS (ICNTL(7)=1, ICNTL(28)=1) gegeben.
//
// Annahmen: FESystem( FE_Q<dim>(k), dim, FE_Q<dim>(1), 1 ), Matrix = TrilinosWrappers::SparseMatrix
// mit zusammenhaengenden Zeilenbereichen pro Rank (wie fuer MumpsReuse), alle DoFs eines Punkts
// gehoeren demselben Rank (in deal.II der Fall; wird geprueft).
//
// Mit gemocktem MPI/deal.II/METIS getestet (1 Rank); mit echtem deal.II/METIS/MPI NICHT kompiliert.
#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/trilinos_index_access.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <metis.h>
#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace NodeND
{
  struct Options
  {
    int  metis_seed         = -1;    // >= 0: METIS_OPTION_SEED
    bool symmetrize         = true;  // A+A^T bilden (billig; METIS verlangt symmetrischen Graphen)
    bool invert_permutation = false; // falls PERM_IN-Konvention von MUMPS umgekehrt sein sollte
    bool verbose            = true;
  };

  // Lokaler Teil des Knotengraphen eines Ranks. Alle Zahlen sind 0-basierte GLOBALE DoF-Indizes;
  // ein Knoten wird durch seinen Repraesentanten (DoF-Index von u_x) identifiziert.
  struct LocalNodeGraph
  {
    int              first_row = 0;  // globaler Index der ersten eigenen Zeile
    std::vector<int> rep_of_row;     // pro eigener Zeile: Repraesentant ihres Knotens
    std::vector<int> node_reps;      // sortierte, eindeutige Repraesentanten der eigenen Knoten
    std::vector<int> node_deg;       // pro eigenem Knoten: Anzahl Nachbar-Repraesentanten
    std::vector<int> adj_reps;       // Nachbar-Repraesentanten, knotenweise hintereinander
  };

  namespace detail
  {
    inline double vm_hwm_gb()
    {
      std::ifstream f("/proc/self/status");
      std::string   line;
      while (std::getline(f, line))
        if (line.rfind("VmHWM:", 0) == 0)
          {
            long kb = 0;
            std::sscanf(line.c_str(), "VmHWM: %ld kB", &kb);
            return kb / 1048576.0;
          }
      return 0.0;
    }

    inline void send_ints(const std::vector<int> &v, int dest, int tag, MPI_Comm comm)
    {
      AssertThrow(v.size() < static_cast<std::size_t>(INT_MAX),
                  dealii::ExcMessage("NodeND: Nachricht zu gross fuer MPI-Zaehler."));
      MPI_Send(v.data(), static_cast<int>(v.size()), MPI_INT, dest, tag, comm);
    }

    inline void recv_ints(std::vector<int> &v, std::size_t n, int src, int tag, MPI_Comm comm)
    {
      AssertThrow(n < static_cast<std::size_t>(INT_MAX),
                  dealii::ExcMessage("NodeND: Nachricht zu gross fuer MPI-Zaehler."));
      v.resize(n);
      MPI_Recv(v.data(), static_cast<int>(n), MPI_INT, src, tag, comm, MPI_STATUS_IGNORE);
    }

    // Selbstkanten entfernen, optional A+A^T, Zeilen sortieren und Duplikate entfernen.
    // Eingabe/Ausgabe: CSR mit 64-bit-Offsets.
    inline void clean_graph(const long long                  N,
                            std::vector<unsigned long long> &xadj,
                            std::vector<int>                &adj,
                            const bool                       symmetrize)
    {
      std::vector<unsigned long long> cnt(N + 1, 0);
      for (long long i = 0; i < N; ++i)
        for (unsigned long long e = xadj[i]; e < xadj[i + 1]; ++e)
          {
            const int j = adj[e];
            if (j == i)
              continue;
            ++cnt[i + 1];
            if (symmetrize)
              ++cnt[j + 1];
          }
      for (long long i = 0; i < N; ++i)
        cnt[i + 1] += cnt[i];

      std::vector<int> adj2(cnt[N]);
      {
        std::vector<unsigned long long> cur(cnt.begin(), cnt.end() - 1);
        for (long long i = 0; i < N; ++i)
          for (unsigned long long e = xadj[i]; e < xadj[i + 1]; ++e)
            {
              const int j = adj[e];
              if (j == i)
                continue;
              adj2[cur[i]++] = j;
              if (symmetrize)
                adj2[cur[j]++] = static_cast<int>(i);
            }
      }
      std::vector<int>().swap(adj); // Eingabe freigeben

      std::vector<int> newlen(N);
#pragma omp parallel for schedule(static)
      for (long long i = 0; i < N; ++i)
        {
          auto first = adj2.begin() + cnt[i];
          auto last  = adj2.begin() + cnt[i + 1];
          std::sort(first, last);
          newlen[i] = static_cast<int>(std::unique(first, last) - first);
        }

      unsigned long long w = 0;
      xadj.assign(N + 1, 0);
      for (long long i = 0; i < N; ++i)
        {
          xadj[i] = w;
          for (int k = 0; k < newlen[i]; ++k)
            adj2[w++] = adj2[cnt[i] + k]; // w <= cnt[i]: Kopie nach vorne ist sicher
        }
      xadj[N] = w;
      adj2.resize(w);
      adj.swap(adj2);
    }
  } // namespace detail


  // Kollektiv. Liefert auf Rank 0 PERM_IN (1-basiert, Laenge n_total), auf den anderen Ranks leer.
  inline std::vector<int> order_node_graph(const LocalNodeGraph &g,
                                           const int             n_total,
                                           MPI_Comm              comm,
                                           const Options        &opt = Options())
  {
    int rank = 0, np = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    const long long head[4] = {g.first_row,
                               static_cast<long long>(g.rep_of_row.size()),
                               static_cast<long long>(g.node_reps.size()),
                               static_cast<long long>(g.adj_reps.size())};
    if (rank != 0)
      {
        MPI_Send(head, 4, MPI_LONG_LONG, 0, 11, comm);
        detail::send_ints(g.rep_of_row, 0, 12, comm);
        detail::send_ints(g.node_reps, 0, 13, comm);
        detail::send_ints(g.node_deg, 0, 14, comm);
        detail::send_ints(g.adj_reps, 0, 15, comm);
        return {};
      }

    // ---------------- Rank 0 ----------------
    const double t0 = MPI_Wtime();

    // Phase A: Zuordnung DoF -> Knoten
    struct Piece
    {
      long long        first_row = 0, n_adj = 0;
      std::vector<int> rep_of_row, node_reps, node_deg;
    };
    std::vector<Piece> pieces(np);
    for (int p = 0; p < np; ++p)
      {
        Piece &pc = pieces[p];
        if (p == 0)
          {
            pc.first_row  = g.first_row;
            pc.n_adj      = head[3];
            pc.rep_of_row = g.rep_of_row;
            pc.node_reps  = g.node_reps;
            pc.node_deg   = g.node_deg;
          }
        else
          {
            long long h[4];
            MPI_Recv(h, 4, MPI_LONG_LONG, p, 11, comm, MPI_STATUS_IGNORE);
            pc.first_row = h[0];
            pc.n_adj     = h[3];
            detail::recv_ints(pc.rep_of_row, h[1], p, 12, comm);
            detail::recv_ints(pc.node_reps, h[2], p, 13, comm);
            detail::recv_ints(pc.node_deg, h[2], p, 14, comm);
          }
      }

    std::vector<long long> node_off(np + 1, 0);
    for (int p = 0; p < np; ++p)
      node_off[p + 1] = node_off[p] + static_cast<long long>(pieces[p].node_reps.size());
    const long long N = node_off[np];
    AssertThrow(N > 0 && N < INT_MAX, dealii::ExcMessage("NodeND: Knotenzahl ungueltig."));

    std::vector<int> node_of_dof(n_total, -1);
    std::vector<int> weight(N, 0);
    for (int p = 0; p < np; ++p)
      for (std::size_t k = 0; k < pieces[p].node_reps.size(); ++k)
        {
          const int rep = pieces[p].node_reps[k];
          AssertThrow(rep >= 0 && rep < n_total && node_of_dof[rep] == -1,
                      dealii::ExcMessage("NodeND: Repraesentant doppelt oder ausserhalb. Die DoFs "
                                         "eines Punkts liegen offenbar auf verschiedenen Ranks."));
          node_of_dof[rep] = static_cast<int>(node_off[p] + k);
        }
    for (int p = 0; p < np; ++p)
      for (std::size_t i = 0; i < pieces[p].rep_of_row.size(); ++i)
        {
          const long long d    = pieces[p].first_row + static_cast<long long>(i);
          const int       rep  = pieces[p].rep_of_row[i];
          AssertThrow(d >= 0 && d < n_total && rep >= 0 && rep < n_total,
                      dealii::ExcMessage("NodeND: Index ausserhalb."));
          const int node = node_of_dof[rep];
          AssertThrow(node >= node_off[p] && node < node_off[p + 1],
                      dealii::ExcMessage("NodeND: Repraesentant gehoert nicht zum eigenen Rank."));
          node_of_dof[d] = node;
          ++weight[node];
        }
    for (int d = 0; d < n_total; ++d)
      AssertThrow(node_of_dof[d] >= 0, dealii::ExcMessage("NodeND: DoF ohne Knoten."));

    // Phase B: Knotengraph einsammeln
    std::vector<unsigned long long> xadj(N + 1, 0);
    for (int p = 0; p < np; ++p)
      for (std::size_t k = 0; k < pieces[p].node_deg.size(); ++k)
        xadj[node_off[p] + k + 1] = pieces[p].node_deg[k];
    for (long long i = 0; i < N; ++i)
      xadj[i + 1] += xadj[i];
    std::vector<int> adj(xadj[N]);
    {
      std::vector<int> buf;
      for (int p = 0; p < np; ++p)
        {
          const std::vector<int> *src = nullptr;
          if (p == 0)
            src = &g.adj_reps;
          else
            {
              detail::recv_ints(buf, pieces[p].n_adj, p, 15, comm);
              src = &buf;
            }
          AssertThrow(xadj[node_off[p + 1]] - xadj[node_off[p]] == src->size(),
                      dealii::ExcMessage("NodeND: Adjazenzlaenge passt nicht zu node_deg."));
          unsigned long long pos = xadj[node_off[p]];
          for (const int r : *src)
            adj[pos++] = node_of_dof[r];
        }
    }
    pieces.clear();
    pieces.shrink_to_fit();
    const double t1 = MPI_Wtime();

    detail::clean_graph(N, xadj, adj, opt.symmetrize);
    const double t2 = MPI_Wtime();
    AssertThrow(xadj[N] < static_cast<unsigned long long>(INT_MAX),
                dealii::ExcMessage("NodeND: Knotengraph hat mehr als 2^31 Eintraege "
                                   "(32-bit-METIS nicht moeglich)."));

    // Phase C: METIS_NodeND mit Vertexgewichten
    std::vector<idx_t> xadj_m(xadj.begin(), xadj.end());
    std::vector<unsigned long long>().swap(xadj);
    std::vector<idx_t> adj_m(adj.begin(), adj.end());
    std::vector<int>().swap(adj);
    std::vector<idx_t> vwgt(weight.begin(), weight.end());

    idx_t nvtxs = static_cast<idx_t>(N);
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    if (opt.metis_seed >= 0)
      options[METIS_OPTION_SEED] = opt.metis_seed;
    std::vector<idx_t> perm(N), iperm(N);
    const int ierr = METIS_NodeND(&nvtxs, xadj_m.data(), adj_m.data(), vwgt.data(), options,
                                  perm.data(), iperm.data());
    AssertThrow(ierr == METIS_OK, dealii::ExcMessage("NodeND: METIS_NodeND fehlgeschlagen."));
    const double t3 = MPI_Wtime();
    const long long n_edges = static_cast<long long>(xadj_m[N]);
    std::vector<idx_t>().swap(xadj_m);
    std::vector<idx_t>().swap(adj_m);

    // Phase D: Expansion auf DoFs. iperm = alt->neu (Position des Knotens in der Eliminationsfolge).
    std::vector<int> order(N, -1); // Knoten an Position q
    for (long long k = 0; k < N; ++k)
      {
        const long long q = iperm[k];
        AssertThrow(q >= 0 && q < N && order[q] == -1,
                    dealii::ExcMessage("NodeND: METIS-iperm ist keine Permutation."));
        order[q] = static_cast<int>(k);
      }
    std::vector<int> start(N);
    long long        acc = 0;
    for (long long q = 0; q < N; ++q)
      {
        start[order[q]] = static_cast<int>(acc);
        acc += weight[order[q]];
      }
    AssertThrow(acc == n_total, dealii::ExcInternalError());

    std::vector<int> perm_in(n_total);
    {
      std::vector<int> fill(N, 0);
      for (int d = 0; d < n_total; ++d)
        {
          const int k = node_of_dof[d];
          perm_in[d]  = start[k] + fill[k]++ + 1; // 1-basiert
        }
    }
    if (opt.invert_permutation)
      {
        std::vector<int> inv(n_total);
        for (int d = 0; d < n_total; ++d)
          inv[perm_in[d] - 1] = d + 1;
        perm_in.swap(inv);
      }
    const double t4 = MPI_Wtime();

    if (opt.verbose)
      {
        int wmin = INT_MAX, wmax = 0;
        for (const int w : weight)
          {
            wmin = std::min(wmin, w);
            wmax = std::max(wmax, w);
          }
        std::cout << "[NODE-ND] n=" << n_total << " Knoten=" << N << " Kanten(gerichtet)="
                  << n_edges << " Gewicht min/max=" << wmin << "/" << wmax << "\n"
                  << "[NODE-ND] Zeit: sammeln " << (t1 - t0) << " s, bereinigen " << (t2 - t1)
                  << " s, METIS " << (t3 - t2) << " s, Expansion " << (t4 - t3)
                  << " s; Rank-0 Speicherspitze " << detail::vm_hwm_gb() << " GB" << std::endl;
      }
    return perm_in;
  }


  // ---------- deal.II-abhaengiger Teil ----------
  template <int dim, int spacedim>
  LocalNodeGraph build_local_node_graph(const dealii::DoFHandler<dim, spacedim>      &dh,
                                        const dealii::TrilinosWrappers::SparseMatrix &A)
  {
    const auto &fe = dh.get_fe();
    AssertThrow(fe.n_base_elements() == 2 && fe.element_multiplicity(0) == dim &&
                  fe.element_multiplicity(1) == 1 && fe.base_element(1).degree == 1,
                dealii::ExcMessage("NodeND: erwartet FESystem(FE_Q(k)^dim, FE_Q(1))."));
    const unsigned int dpc     = fe.n_dofs_per_cell();
    const unsigned int n_unode = fe.base_element(0).n_dofs_per_cell();
    AssertThrow(fe.base_element(1).n_dofs_per_cell() <= n_unode,
                dealii::ExcMessage("NodeND: p-Element hat mehr Knoten als u-Element."));

    // 1) (DoF, Repraesentant) fuer alle DoFs aller nicht-kuenstlichen Zellen.
    //    Lokaler u-Knoten j <-> p-Knoten j (FE_Q: Ecken zuerst, gleiche Reihenfolge).
    std::vector<std::pair<int, int>>           pairs;
    std::vector<dealii::types::global_dof_index> ldi(dpc);
    std::vector<int>                           ux(n_unode);
    for (const auto &cell : dh.active_cell_iterators())
      {
        if (cell->is_artificial())
          continue;
        cell->get_dof_indices(ldi);
        std::fill(ux.begin(), ux.end(), -1);
        for (unsigned int i = 0; i < dpc; ++i)
          {
            const auto b = fe.system_to_base_index(i);
            if (b.first.first == 0 && b.first.second == 0)
              ux[b.second] = static_cast<int>(ldi[i]);
          }
        for (unsigned int i = 0; i < dpc; ++i)
          {
            const auto b = fe.system_to_base_index(i);
            AssertThrow(ux[b.second] >= 0, dealii::ExcInternalError());
            pairs.emplace_back(static_cast<int>(ldi[i]), ux[b.second]);
          }
      }
    std::sort(pairs.begin(), pairs.end());
    {
      std::size_t w = 0;
      for (std::size_t i = 0; i < pairs.size(); ++i)
        {
          if (w > 0 && pairs[w - 1].first == pairs[i].first)
            {
              AssertThrow(pairs[w - 1].second == pairs[i].second,
                          dealii::ExcMessage("NodeND: uneindeutiger Repraesentant."));
              continue;
            }
          pairs[w++] = pairs[i];
        }
      pairs.resize(w);
    }
    auto rep_of = [&](const long long d) -> int {
      auto it = std::lower_bound(pairs.begin(),
                                 pairs.end(),
                                 std::make_pair(static_cast<int>(d), INT_MIN));
      AssertThrow(it != pairs.end() && it->first == d,
                  dealii::ExcMessage("NodeND: DoF der Matrix ist in keiner lokalen/Ghost-Zelle."));
      return it->second;
    };

    // 2) Matrix: Repraesentanten der Spalten und Zeilen
    const auto &M      = A.trilinos_matrix();
    const int   n_rows = M.NumMyRows();
    const int   n_cols = M.NumMyCols();

    LocalNodeGraph g;
    g.first_row = n_rows > 0 ? static_cast<int>(dealii::TrilinosWrappers::global_row_index(M, 0)) : 0;
    g.rep_of_row.resize(n_rows);
    for (int r = 0; r < n_rows; ++r)
      {
        const long long gr = dealii::TrilinosWrappers::global_row_index(M, r);
        AssertThrow(gr == static_cast<long long>(g.first_row) + r,
                    dealii::ExcMessage("NodeND: eigene Zeilen nicht zusammenhaengend."));
        g.rep_of_row[r] = rep_of(gr);
      }
    std::vector<int> rep_of_col(n_cols);
    for (int c = 0; c < n_cols; ++c)
      rep_of_col[c] = rep_of(dealii::TrilinosWrappers::global_column_index(M, c));

    // 3) Knoten = eindeutige Repraesentanten der eigenen Zeilen
    g.node_reps = g.rep_of_row;
    std::sort(g.node_reps.begin(), g.node_reps.end());
    g.node_reps.erase(std::unique(g.node_reps.begin(), g.node_reps.end()), g.node_reps.end());
    const int n_nodes = static_cast<int>(g.node_reps.size());

    // Zeilen nach Knoten gruppieren (Counting Sort)
    std::vector<int> node_of_row(n_rows), rows_ptr(n_nodes + 1, 0), rows_sorted(n_rows);
    for (int r = 0; r < n_rows; ++r)
      {
        node_of_row[r] = static_cast<int>(
          std::lower_bound(g.node_reps.begin(), g.node_reps.end(), g.rep_of_row[r]) -
          g.node_reps.begin());
        ++rows_ptr[node_of_row[r] + 1];
      }
    for (int k = 0; k < n_nodes; ++k)
      rows_ptr[k + 1] += rows_ptr[k];
    {
      std::vector<int> cur(rows_ptr.begin(), rows_ptr.end() - 1);
      for (int r = 0; r < n_rows; ++r)
        rows_sorted[cur[node_of_row[r]]++] = r;
    }

    // 4) Nachbar-Repraesentanten je Knoten (Vereinigung ueber alle Zeilen des Knotens)
    g.node_deg.assign(n_nodes, 0);
    std::vector<int> tmp;
    for (int k = 0; k < n_nodes; ++k)
      {
        tmp.clear();
        for (int s = rows_ptr[k]; s < rows_ptr[k + 1]; ++s)
          {
            int     num_entries;
            double *values;
            int    *local_cols;
            const int ierr = M.ExtractMyRowView(rows_sorted[s], num_entries, values, local_cols);
            AssertThrow(ierr == 0, dealii::ExcInternalError());
            for (int j = 0; j < num_entries; ++j)
              tmp.push_back(rep_of_col[local_cols[j]]);
          }
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        for (const int r : tmp)
          if (r != g.node_reps[k])
            {
              g.adj_reps.push_back(r);
              ++g.node_deg[k];
            }
      }
    return g;
  }


  // Bequemer Einstieg: kollektiv; PERM_IN (1-basiert) nur auf Rank 0.
  template <int dim, int spacedim>
  std::vector<int> compute_permutation(const dealii::DoFHandler<dim, spacedim>      &dh,
                                       const dealii::TrilinosWrappers::SparseMatrix &A,
                                       MPI_Comm                                      comm,
                                       const Options                                &opt = Options())
  {
    AssertThrow(A.m() < static_cast<std::size_t>(INT_MAX),
                dealii::ExcMessage("NodeND: n >= 2^31 nicht unterstuetzt."));
    const LocalNodeGraph g = build_local_node_graph(dh, A);
    return order_node_graph(g, static_cast<int>(A.m()), comm, opt);
  }
} // namespace NodeND
