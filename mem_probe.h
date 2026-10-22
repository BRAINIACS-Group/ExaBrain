// mem_probe.h
// Misst den tatsaechlich belegten Hauptspeicher (VmRSS aus /proc/self/status) an beliebigen Stellen,
// aufsummiert pro Knoten (MPI_Comm_split_type SHARED) und ueber alle Knoten. Das ist dieselbe Groesse
// wie "mem_used" im ClusterCockpit (residenter Speicher aller Prozesse eines Knotens).
// Kollektiv aufrufen (alle Ranks):   MemProbe::report("nach system_setup", mpi_communicator);
// Ausgabe (Rank 0):  [MEM] <tag>  Knoten RSS min/avg/max = a/b/c GiB  Summe = s GiB  Rank0 = r GiB  HWM(Knoten, max) = h GiB
#pragma once

#include <mpi.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace MemProbe
{
  inline void read_status(double &rss_gib, double &hwm_gib)
  {
    rss_gib = hwm_gib = 0.0;
    std::ifstream f("/proc/self/status");
    std::string   line;
    while (std::getline(f, line))
      {
        long kb = 0;
        if (line.rfind("VmRSS:", 0) == 0 && std::sscanf(line.c_str(), "VmRSS: %ld kB", &kb) == 1)
          rss_gib = kb / 1048576.0;
        else if (line.rfind("VmHWM:", 0) == 0 && std::sscanf(line.c_str(), "VmHWM: %ld kB", &kb) == 1)
          hwm_gib = kb / 1048576.0;
      }
  }

  inline void report(const char *tag, MPI_Comm comm)
  {
    MPI_Barrier(comm);
    double rss, hwm;
    read_status(rss, hwm);

    MPI_Comm nc;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nc);
    const double loc[2] = {rss, hwm};
    double       node[2];
    MPI_Allreduce(loc, node, 2, MPI_DOUBLE, MPI_SUM, nc); // Summe ueber die Ranks eines Knotens
    int nrank;
    MPI_Comm_rank(nc, &nrank);
    MPI_Comm_free(&nc);

    const bool   lead = (nrank == 0); // ein Vertreter pro Knoten
    double       mn   = lead ? node[0] : 1e300;
    double       mx   = lead ? node[0] : 0.0;
    double       hmx  = lead ? node[1] : 0.0;
    const double sums[2] = {lead ? node[0] : 0.0, lead ? 1.0 : 0.0}; // RSS-Summe, Knotenzahl
    double       s[2];
    MPI_Allreduce(sums, s, 2, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &mn, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(MPI_IN_PLACE, &hmx, 1, MPI_DOUBLE, MPI_MAX, comm);

    int rank;
    MPI_Comm_rank(comm, &rank);
    double r0 = rss;
    MPI_Bcast(&r0, 1, MPI_DOUBLE, 0, comm);
    if (rank == 0)
      std::cout << "[MEM] " << tag << "  Knoten RSS min/avg/max = " << mn << "/" << s[0] / s[1] << "/"
                << mx << " GiB  Summe = " << s[0] << " GiB (" << static_cast<int>(s[1])
                << " Knoten)  Rank0 = " << r0 << " GiB  HWM(Knoten, max) = " << hmx << " GiB"
                << std::endl;
  }
} // namespace MemProbe
