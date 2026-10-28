// diag_gap.h
#include <mpi.h>
#include <iostream>
#include <mutex>
#include <vector>

namespace DiagGap
{
  struct alignas(64) Local { unsigned long b[8] = {0,0,0,0,0,0,0,0}; };

  inline std::mutex           reg_mutex;
  inline std::vector<Local *> registry;   // je ein Zaehlerblock pro Thread

  inline Local &local()
  {
    // wird pro Thread genau einmal angelegt und eingetragen (nie freigegeben, 64 Byte)
    thread_local Local *p = []() {
      Local *l = new Local;
      std::lock_guard<std::mutex> g(reg_mutex);
      registry.push_back(l);
      return l;
    }();
    return *p;
  }

  // Klassen: <1e-14, <1e-12, <1e-10, <1e-8, <1e-6, <1e-4, <1e-2, sonst
  inline void count(const double gap)
  {
    static const double edge[7] = {1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
    unsigned k = 0;
    while (k < 7 && gap >= edge[k]) ++k;
    ++local().b[k];
  }

  // NUR aufrufen, wenn keine Assemblierungs-Threads mehr laufen (nach WorkStream::run)
  inline void report_and_reset(MPI_Comm comm, const char *tag)
  {
    unsigned long loc[8] = {0,0,0,0,0,0,0,0}, glob[8];
    {
      std::lock_guard<std::mutex> g(reg_mutex);
      for (Local *l : registry)
        for (int k = 0; k < 8; ++k) { loc[k] += l->b[k]; l->b[k] = 0; }
    }
    MPI_Allreduce(loc, glob, 8, MPI_UNSIGNED_LONG, MPI_SUM, comm);
    int r; MPI_Comm_rank(comm, &r);
    if (r == 0)
      {
        std::cout << "[GAP] " << tag;
        for (int k = 0; k < 8; ++k) std::cout << " " << glob[k];
        std::cout << std::endl;
      }
  }
}
