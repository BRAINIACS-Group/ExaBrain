// mumps_reuse.h
// Eigener MUMPS-Wrapper fuer dealii::TrilinosWrappers::SparseMatrix (Epetra, unsymmetrisch).
// Gegenueber deal.II SparseDirectMUMPS:
//   * Analyse (JOB=1) nur einmal, danach pro initialize() nur JOB=2 (Faktorisierung)
//   * optional mehrfache Analyse beim ersten Aufruf (Auswahl nach RINFOG(1)), fuer ParMETIS
//   * Tripel irn/jcn werden nur einmal aufgebaut, danach nur die Werte aktualisiert
//   * Loesung wird per MPI_Scatterv verteilt (Zaehler/Offsets einmal berechnet)
//   * Zeitmessung von MUMPS-Solve und Umverteilung
// NICHT kompiliert/getestet. Benoetigt dmumps_c.h im Include-Pfad (wie dein Projekt schon).
#pragma once

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/lac/trilinos_index_access.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include "mem_probe.h"

#include <dmumps_c.h>
#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <utility>
#include <string>
#include <vector>

class MumpsReuse
{
public:
  // Gleiche Feldnamen wie dealii::SparseDirectMUMPS::AdditionalData, damit
  // data_mumps.blr.lowrank_threshold usw. unveraendert bleiben.
  struct AdditionalData
  {
    struct BLR
    {
      double lowrank_threshold = 1e-4;
      bool   blr_ucfs          = true;
    };

    bool blr_factorization = false;
    bool output_details    = false;
    BLR  blr;
  };

  MumpsReuse(const AdditionalData &data, const MPI_Comm &comm)
    : additional_data(data)
    , comm(comm)
  {
    id.job          = -1;
    id.par          = 1;
    id.sym          = 0;
    id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(comm);
    dmumps_c(&id);

    if (!additional_data.output_details)
      {
        id.icntl[0] = -1;
        id.icntl[1] = -1;
        id.icntl[2] = -1;
        id.icntl[3] = 0;
      }
    if (additional_data.blr_factorization)
      {
        id.icntl[34] = 2;
        if (additional_data.blr.blr_ucfs)
          id.icntl[35] = 1;
        AssertThrow(additional_data.blr.lowrank_threshold > 0,
                    dealii::ExcMessage("Lowrank threshold must be positive."));
        id.cntl[6] = additional_data.blr.lowrank_threshold;
      }
  }

  ~MumpsReuse()
  {
    id.job = -2;
    dmumps_c(&id);
  }

  MumpsReuse(const MumpsReuse &)            = delete;
  MumpsReuse &operator=(const MumpsReuse &) = delete;

  // Zugriff wie beim bisherigen Wrapper: IC(i) = icntl[i-1], RC(i) = cntl[i-1]
  int    *get_icntl() { return id.icntl; }
  double *get_cntl() { return id.cntl; }

  // tries: Anzahl der Analysen beim ersten initialize().
  // accept_rinfog1 > 0: Abbruch, sobald RINFOG(1) <= Schwelle. 0: immer alle 'tries' Analysen.
  // Es bleibt die LETZTE Analyse aktiv (MUMPS speichert keine beste).
  void set_analysis_policy(const unsigned int tries, const double accept_rinfog1 = 0.0)
  {
    max_tries     = std::max(1u, tries);
    accept_flops  = accept_rinfog1;
  }

  // Nach jeder Gitteraenderung / neuem system_setup aufrufen.
  void reset_analysis()
  {
    if (use_user_ordering)
      {
        id.icntl[6]  = saved_icntl7;  // vom Nutzer gesetzte Werte wiederherstellen
        id.icntl[27] = saved_icntl28;
      }
    use_user_ordering = false;
    perm_in.clear();
    perm_in.shrink_to_fit();
    analysis_done = false;
  }

  bool analysis_is_done() const { return analysis_done; }

  // Speichermessung (MemProbe) nach Triplet-Aufbau, Analyse und Faktorisierung ein/aus.
  void set_memory_probe(const bool on) { probe_memory = on; }

  // Eigene Pivotreihenfolge (z.B. NodeND::compute_permutation). Nur auf Rank 0 nicht leer,
  // 1-basiert, Laenge n. Vor dem ersten initialize() nach reset_analysis() aufrufen.
  // Setzt bei der Analyse ICNTL(7)=1 und ICNTL(28)=1 (sequentielle Analyse).
  void set_user_ordering(std::vector<MUMPS_INT> &&perm)
  {
    perm_in           = std::move(perm);
    use_user_ordering = true;
  }

  void initialize(const dealii::TrilinosWrappers::SparseMatrix &matrix)
  {
    if (!analysis_done)
      {
        probe("vor Triplet-Aufbau");
        build_pattern(matrix);
        probe("nach Triplet-Aufbau");
        run_analysis();
        probe("nach Analyse");
        analysis_done = true;
      }
    else
      update_values(matrix);

    id.job = 2; // nur Faktorisierung
    dmumps_c(&id);
    check_status("factorization");
    probe("nach Faktorisierung");
  }

  // FGMRES-Praekonditionierer-Schnittstelle
  void vmult(dealii::TrilinosWrappers::MPI::Vector       &dst,
             const dealii::TrilinosWrappers::MPI::Vector &src) const
  {
    AssertThrow(analysis_done, dealii::ExcMessage("MumpsReuse: initialize() fehlt."));
    const bool root = (dealii::Utilities::MPI::this_mpi_process(comm) == 0);

    id.rhs_loc = const_cast<double *>(src.begin());
    if (root)
      {
        if (rhs.size() != static_cast<std::size_t>(n))
          rhs.resize(n);
        id.rhs = rhs.data();
      }

    const double t0 = MPI_Wtime();
    id.job          = 3;
    dmumps_c(&id);
    const double t1 = MPI_Wtime();

    MPI_Scatterv(root ? rhs.data() : nullptr,
                 counts.data(),
                 displs.data(),
                 MPI_DOUBLE,
                 dst.begin(),
                 local_n,
                 MPI_DOUBLE,
                 0,
                 comm);
    const double t2 = MPI_Wtime();

    t_mumps_solve += t1 - t0;
    t_redistribute += t2 - t1;
    ++n_solves;
  }

  // Am Ende von solve_linear_system oder des Programms aufrufen (kollektiv).
  void print_solve_stats(dealii::ConditionalOStream &pcout) const
  {
    const auto a = dealii::Utilities::MPI::min_max_avg(t_mumps_solve, comm);
    const auto b = dealii::Utilities::MPI::min_max_avg(t_redistribute, comm);
    pcout << "[MUMPS-SOLVE] n=" << n_solves << "  JOB=3 avg/max " << a.avg << "/" << a.max
          << " s  redistribute avg/max " << b.avg << "/" << b.max << " s" << std::endl;
  }

private:
  void build_pattern(const dealii::TrilinosWrappers::SparseMatrix &matrix)
  {
    AssertThrow(matrix.n() == matrix.m(), dealii::ExcMessage("Matrix must be square."));
    int result;
    MPI_Comm_compare(comm, matrix.get_mpi_communicator(), &result);
    AssertThrow(result == MPI_IDENT, dealii::ExcMessage("Communicator mismatch."));

    n    = matrix.n();
    id.n = n;

    const auto &M       = matrix.trilinos_matrix();
    const int   n_rows  = M.NumMyRows();
    const auto  n_local = static_cast<std::size_t>(M.NumMyNonzeros());

    irn.assign(n_local, 0);
    jcn.assign(n_local, 0);
    a.assign(n_local, 0.0);
    irhs_loc.assign(n_rows, 0);

    std::size_t k = 0;
    for (int r = 0; r < n_rows; ++r)
      {
        int     num_entries;
        double *values;
        int    *local_cols;
        const int ierr = M.ExtractMyRowView(r, num_entries, values, local_cols);
        AssertThrow(ierr == 0, dealii::ExcInternalError());
        const auto grow = dealii::TrilinosWrappers::global_row_index(M, r);
        for (int j = 0; j < num_entries; ++j)
          {
            irn[k] = grow + 1;
            jcn[k] = dealii::TrilinosWrappers::global_column_index(M, local_cols[j]) + 1;
            a[k]   = values[j];
            ++k;
          }
        irhs_loc[r] = grow + 1;
      }
    AssertThrow(k == n_local, dealii::ExcInternalError());

    // Zeilen muessen pro Rank zusammenhaengend aufsteigend sein (wie im deal.II-Wrapper).
    for (int r = 1; r < n_rows; ++r)
      AssertThrow(irhs_loc[r] == irhs_loc[r - 1] + 1,
                  dealii::ExcMessage("Locally owned rows are not contiguous."));

    id.icntl[17] = 3;   // verteilte Matrix
    id.nnz       = matrix.n_nonzero_elements();
    id.nnz_loc   = n_local;
    id.irn_loc   = irn.data();
    id.jcn_loc   = jcn.data();
    id.a_loc     = a.data();
    id.irhs_loc  = irhs_loc.data();

    id.icntl[19] = 10;  // verteilte rhs
    id.icntl[20] = 0;   // zentrale Loesung auf Rank 0
    id.nrhs      = 1;
    id.lrhs_loc  = n;
    id.nloc_rhs  = n_rows;

    // Zaehler/Offsets fuer Scatterv, einmalig
    local_n = n_rows;
    const int np = dealii::Utilities::MPI::n_mpi_processes(comm);
    counts.assign(np, 0);
    displs.assign(np, 0);
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
    for (int p = 1; p < np; ++p)
      displs[p] = displs[p - 1] + counts[p - 1];
    n_nnz_check = n_local;
  }

  void update_values(const dealii::TrilinosWrappers::SparseMatrix &matrix)
  {
    const auto &M = matrix.trilinos_matrix();
    AssertThrow(static_cast<std::size_t>(M.NumMyNonzeros()) == n_nnz_check,
                dealii::ExcMessage("Sparsity pattern changed since the analysis. "
                                   "Call reset_analysis()."));
    std::size_t k = 0;
    for (int r = 0; r < M.NumMyRows(); ++r)
      {
        int     num_entries;
        double *values;
        int    *local_cols;
        const int ierr = M.ExtractMyRowView(r, num_entries, values, local_cols);
        AssertThrow(ierr == 0, dealii::ExcInternalError());
        for (int j = 0; j < num_entries; ++j)
          a[k++] = values[j];
      }
    AssertThrow(k == n_nnz_check, dealii::ExcInternalError());
    id.a_loc = a.data(); // Zeiger unveraendert, zur Sicherheit
  }

  void run_analysis()
  {
    const bool root = (dealii::Utilities::MPI::this_mpi_process(comm) == 0);
    if (use_user_ordering)
      {
        saved_icntl7  = id.icntl[6];
        saved_icntl28 = id.icntl[27];
        id.icntl[6]   = 1;
        id.icntl[27]  = 1;
        if (root)
          {
            AssertThrow(perm_in.size() == static_cast<std::size_t>(n),
                        dealii::ExcMessage("MumpsReuse: PERM_IN hat falsche Laenge."));
            id.perm_in = perm_in.data();
          }
      }
    for (unsigned int attempt = 0; attempt < max_tries; ++attempt)
      {
        id.job = 1; // nur Analyse
        dmumps_c(&id);
        check_status("analysis");

        // RINFOG(1), INFOG(5), INFOG(7) (verwendete Ordering-Methode, 1 = PERM_IN)
        double buf[3] = {id.rinfog[0], static_cast<double>(id.infog[4]),
                         static_cast<double>(id.infog[6])};
        MPI_Bcast(buf, 3, MPI_DOUBLE, 0, comm);
        if (root)
          std::cout << "[MUMPS-ANALYSIS] try " << attempt << ": RINFOG(1)=" << buf[0]
                    << " INFOG(5)=" << buf[1] << " INFOG(7)=" << buf[2] << std::endl;

        if (accept_flops > 0.0 && buf[0] <= accept_flops)
          break;
        if (attempt + 1 == max_tries && accept_flops > 0.0 && buf[0] > accept_flops && root)
          std::cout << "[MUMPS-ANALYSIS] WARNING: Schwelle " << accept_flops
                    << " nicht erreicht, letzte Analyse bleibt aktiv." << std::endl;
      }
  }

  void check_status(const char *what) const
  {
    int status = id.infog[0];
    MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_MIN, comm);
    AssertThrow(status >= 0,
                dealii::ExcMessage(std::string("MUMPS ") + what +
                                   " failed, INFOG(1) = " + std::to_string(status) +
                                   (status == -9 ? " (mehr Platz: ICNTL(14) erhoehen)" : "")));
  }

  AdditionalData additional_data;
  MPI_Comm       comm;
  mutable DMUMPS_STRUC_C id;

  void probe(const char *tag) const
  {
    if (probe_memory)
      MemProbe::report(tag, comm);
  }

  bool                   probe_memory      = false;
  bool                   use_user_ordering = false;
  std::vector<MUMPS_INT> perm_in;
  int                    saved_icntl7 = 0, saved_icntl28 = 0;

  bool         analysis_done = false;
  unsigned int max_tries     = 1;
  double       accept_flops  = 0.0;

  MUMPS_INT              n = 0;
  std::vector<MUMPS_INT> irn, jcn, irhs_loc;
  std::vector<double>    a;
  std::size_t            n_nnz_check = 0;

  int              local_n = 0;
  std::vector<int> counts, displs;
  mutable std::vector<double> rhs;

  mutable double       t_mumps_solve  = 0.0;
  mutable double       t_redistribute = 0.0;
  mutable unsigned int n_solves       = 0;
};
