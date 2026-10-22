#ifndef PHASE_PROBE_H
#define PHASE_PROBE_H

// PhaseProbe: Zeitmessung der Abschnitte in assemble_system_one_cell (pro Thread, ohne Locks).
//
// Einschalten zur Laufzeit:   export PHASE_PROBE=1   (und  -x PHASE_PROBE  bei mpirun)
// Ausgeschaltet kostet ein Lap-Aufruf nur einen Vergleich mit einer static-Variable.
//
// Benutzung:
//   PhaseProbe::Lap lap;                    // am Anfang der Zellfunktion
//   ... Code ...
//   lap(PhaseProbe::QP_SOLUTION);           // addiert die Zeit seit dem letzten lap() zur Phase
//   ...
//   lap.count_cell();                       // einmal pro Zelle
//   Am Ende des Laufs (nach der letzten Assembly, ausserhalb jeder parallelen Region):
//   PhaseProbe::report(mpi_communicator);

#include <mpi.h>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace PhaseProbe
{
  enum Phase : int
  {
    RESET_REINIT = 0,   // data.reset, scratch.reset, FEValues::reinit, get_dof_indices
    AD_SEED,            // local_dof_values + diff()
    QP_SOLUTION,        // Akkumulation von Gradienten/Druck ueber alle DoFs
    PRECALC,            // symm_grad_Nx, grad_Nx_p_fluid
    KINEMATICS,         // F, det F, F^-1, p am QP
    MAT_UPDATE,         // update_internal_equilibrium (lokales Newton)
    MAT_STRESS,         // get_tau_E, tau_fluid_vol, converged det F, body force
    DARCY,              // Aliase, grad_p, Sickergeschwindigkeit
    RESIDUAL,           // Residuenschleife ueber die DoFs
    NEUMANN,            // Randterme
    LINEARIZE,          // R.val()/fastAccessDx in cell_matrix
    COPY_GLOBAL,        // optional: Copier (in copy_local_to_global verwenden)
    N_PHASES
  };

  inline const char *name(const int p)
  {
    static const char *n[N_PHASES] = {"reset/reinit", "AD-Seeding", "QP-Loesung (Akkum.)",
                                      "Vorberechnung N", "Kinematik F/Finv", "Material update",
                                      "Material Spannung", "Darcy", "Residuenschleife",
                                      "Neumann", "Linearisierung", "Copier"};
    return n[p];
  }

  constexpr int MAX_THREADS = 512;

  constexpr int N_HIST = 17;   // Iterationen 0..15, Bucket 16 = mehr

  struct alignas(64) Row
  {
    double t[N_PHASES + 1];   // letzte Stelle: Zellzaehler
    double hist[N_HIST];      // Histogramm lokale Newton-Iterationen (Material)
    Row()
    {
      for (double &x : t) x = 0.0;
      for (double &x : hist) x = 0.0;
    }
  };

  inline Row *rows()
  {
    static Row r[MAX_THREADS];
    return r;
  }

  inline bool enabled()
  {
    static const bool e = (std::getenv("PHASE_PROBE") != nullptr);
    return e;
  }

  class Lap
  {
  public:
    Lap() : on(enabled()), row(nullptr)
    {
      if (on)
        {
          row  = &rows()[omp_get_thread_num() % MAX_THREADS];
          last = clock::now();
        }
    }
    void operator()(const int phase)
    {
      if (!on) return;
      const clock::time_point now = clock::now();
      row->t[phase] += std::chrono::duration<double>(now - last).count();
      last = now;
    }
    void count_cell() { if (on) row->t[N_PHASES] += 1.0; }
  private:
    typedef std::chrono::steady_clock clock;
    bool on;
    Row *row;
    clock::time_point last;
  };

  // Im Material nach der lokalen Newton-Schleife aufrufen: PhaseProbe::count_newton(iteration);
  inline void count_newton(const int iterations)
  {
    if (!enabled()) return;
    Row &r = rows()[omp_get_thread_num() % MAX_THREADS];
    r.hist[iterations < 0 ? 0 : (iterations >= N_HIST ? N_HIST - 1 : iterations)] += 1.0;
  }

  inline void reset()
  {
    for (int t = 0; t < MAX_THREADS; ++t)
      {
        for (double &x : rows()[t].t) x = 0.0;
        for (double &x : rows()[t].hist) x = 0.0;
      }
  }

  // Nur aufrufen, wenn keine parallele Region mehr laeuft.
  inline void report(MPI_Comm comm, const char *tag = "")
  {
    if (!enabled()) return;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double loc[N_PHASES + 1];
    for (double &x : loc) x = 0.0;
    int threads_used = 0;
    for (int t = 0; t < MAX_THREADS; ++t)
      {
        bool used = false;
        for (int p = 0; p <= N_PHASES; ++p)
          {
            loc[p] += rows()[t].t[p];
            if (rows()[t].t[p] != 0.0) used = true;
          }
        if (used) ++threads_used;
      }

    double hloc[N_HIST];
    for (double &x : hloc) x = 0.0;
    for (int t = 0; t < MAX_THREADS; ++t)
      for (int b = 0; b < N_HIST; ++b)
        hloc[b] += rows()[t].hist[b];
    double hsum[N_HIST];
    MPI_Allreduce(hloc, hsum, N_HIST, MPI_DOUBLE, MPI_SUM, comm);

    double sum[N_PHASES + 1], mx[N_PHASES + 1];
    MPI_Allreduce(loc, sum, N_PHASES + 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(loc, mx,  N_PHASES + 1, MPI_DOUBLE, MPI_MAX, comm);
    double thr_sum = threads_used, thr_max = threads_used;
    double tu = threads_used;
    MPI_Allreduce(&tu, &thr_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&tu, &thr_max, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (rank != 0) return;

    double total = 0.0;
    for (int p = 0; p < N_PHASES; ++p) total += sum[p];
    const double cells = sum[N_PHASES];

    std::printf("[PHASE] %s  Zellen=%.0f  Raenge=%d  Threads/Rang(max)=%.0f  Summe Thread-Sekunden=%.1f  mittlere Zeit je Thread=%.1f s\n",
                tag, cells, size, thr_max, total, thr_sum > 0 ? total / thr_sum : 0.0);
    std::printf("[PHASE] %-22s %12s %8s %14s %10s\n", "Phase", "Thread-s", "Anteil", "us/Zelle", "Max/Mittel");
    for (int p = 0; p < N_PHASES; ++p)
      {
        const double avg_rank = sum[p] / size;
        std::printf("[PHASE] %-22s %12.2f %7.1f%% %14.2f %10.2f\n", name(p), sum[p],
                    total > 0 ? 100.0 * sum[p] / total : 0.0,
                    cells > 0 ? 1e6 * sum[p] / cells : 0.0,
                    avg_rank > 0 ? mx[p] / avg_rank : 0.0);
      }
    double calls = 0.0, iters = 0.0;
    for (int b = 0; b < N_HIST; ++b) { calls += hsum[b]; iters += hsum[b] * b; }
    if (calls > 0)
      {
        std::printf("[PHASE] lokales Newton: Aufrufe=%.0f  mittlere Iterationen=%.2f  Verteilung (Iter: Anteil):", calls, iters / calls);
        for (int b = 0; b < N_HIST; ++b)
          if (hsum[b] > 0)
            std::printf("  %d%s: %.1f%%", b, b == N_HIST - 1 ? "+" : "", 100.0 * hsum[b] / calls);
        std::printf("\n");
      }
    std::fflush(stdout);
  }
}

#endif
