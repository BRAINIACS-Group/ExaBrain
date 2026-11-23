# Nonlinear poro-viscoelasticity (deal.II)

Finite element code for finite-strain, biphasic (solid + fluid) poro-viscoelastic
brain tissue, based on the deal.II code gallery example by E. Comellas and
J.-P. Pelteret. Everything is in a single source file,
`nonlinear-poro-viscoelasticity.cc`, plus a few helper headers. Each run is
controlled by one parameter file (`.prm`).

This README covers building on macOS (Apple Silicon) and running the two
quarter-cylinder test cases in `parameter_files/`.

---

## 1. Requirements

### deal.II

The code needs a recent **deal.II master (9.8.0-pre)**; deal.II 9.7.1 does not
compile it (`get_cntl` error). Tested with master commit `95600ced65`
(2025-10-23).

deal.II must be configured with:

| Feature | Why |
|---|---|
| `DEAL_II_WITH_MPI` | distributed runs |
| `DEAL_II_WITH_TRILINOS` with Epetra, EpetraExt, **Sacado** | matrices/vectors, automatic differentiation of the material models (checked by `CMakeLists.txt`) |
| `DEAL_II_WITH_MUMPS` | linear solvers. The code also includes `dmumps_c.h` and calls MUMPS **directly**, so deal.II, Trilinos and the code must all use the *same* MUMPS |
| `DEAL_II_WITH_METIS` | used directly (`metis.h`) |
| `DEAL_II_WITH_P4EST` | distributed triangulation |
| `DEAL_II_WITH_BOOST` (with `filesystem`) | used directly |
| `DEAL_II_WITH_LAPACK`, `DEAL_II_WITH_SCALAPACK` | needed by MUMPS |

Also enabled in the tested build (not strictly needed by this code): Kokkos,
TBB, zlib, Taskflow, muparser, UMFPACK. PETSc and SLEPc are **off**.

### MUMPS

Tested: **MUMPS 5.9.1** (`+double +metis +parmetis +scotch +ptscotch +mpi +shared ~int64 +openmp`).
Build MUMPS **with OpenMP**. The `BLRNewt` solver uses MUMPS options
(ICNTL(48)) that fail with a MUMPS built without OpenMP. `MumpsDirect` works
either way.

### OpenMP on macOS

Apple clang has no OpenMP; use Homebrew `libomp`. The code includes `omp.h`, so
the project is compiled with `-Xpreprocessor -fopenmp` and linked with `-lomp`
(see section 3).

---

## 2. Building the dependencies on macOS (tested setup)

Homebrew packages:

```bash
brew install cmake open-mpi openblas libomp gcc
```

(`gcc` provides `gfortran` for MUMPS/ScaLAPACK.)

The rest comes from a [Spack](https://spack.io) environment. Tested `spack.yaml`:

```yaml
spack:
  specs:
  - mumps@5.9.1 +float+double+complex +metis+parmetis+scotch+ptscotch +mpi +shared ~int64 +openmp ~blr_mt fflags="-fopenmp" cflags="-I/opt/homebrew/opt/libomp/include"
  - trilinos@16.2.2 +amesos+amesos2+anasazi+aztec+belos+epetra+epetraext+ifpack+ifpack2+kokkos+ml+muelu+mumps+rol+sacado+superlu-dist+suite-sparse+tpetra+mesquite+explicit_template_instantiation+fortran+mpi+shared +openmp~hdf5~boost~zoltan~zoltan2~stratimikos~thyra~nox~teko~shylu cxxstd=17 build_type=Release ldflags="-L/opt/homebrew/opt/libomp/lib" ldlibs="-lomp"
  - p4est@2.8.7 +mpi ~openmp
  - metis@5.1.0 ~int64
  - parmetis@4.0.3 ~int64
  - netlib-scalapack@2.2.2
  - boost@1.89.0 +atomic+chrono+container+date_time+exception+filesystem+graph+icu+iostreams+locale+log+math+program_options+random+regex+serialization+system+test+thread+timer+wave+multithreaded+shared cxxstd=17
  - intel-tbb
  - zlib-ng
  view: <env-dir>/view
  concretizer:
    unify: true
    reuse: true
  packages:
    all:
      providers:
        mpi: [openmpi]
        blas: [openblas]
        lapack: [openblas]
        scalapack: [netlib-scalapack]
    openmpi:            # use Homebrew's Open MPI
      externals:
      - spec: openmpi@5.0.9
        prefix: /opt/homebrew/opt/open-mpi
      buildable: false
    openblas:           # use Homebrew's OpenBLAS
      externals:
      - spec: openblas@0.3.33
        prefix: /opt/homebrew/opt/openblas
      buildable: false
    llvm-openmp:        # use Homebrew's libomp
      externals:
      - spec: llvm-openmp@22.1.4
        prefix: /opt/homebrew/opt/libomp
      buildable: false
    kokkos:
      require: "+serial +openmp"
```

Adjust the versions in the `externals` entries to what Homebrew actually
installed (`brew list --versions open-mpi openblas libomp`). Then:

```bash
spack env create dealii-deps spack.yaml
spack env activate dealii-deps
spack install
```

**Known pitfall (only for more than 1 thread per MPI rank):** `gfortran` links
MUMPS against GCC's `libgomp`, while the rest of the stack uses LLVM `libomp`.
Mixing both breaks MUMPS with several threads per rank. With 1 thread per rank
(recommended, see section 4) this does not matter. If you need threads, relink
the MUMPS dylibs to `libomp` (`install_name_tool -change .../libgomp.1.dylib
/opt/homebrew/opt/libomp/lib/libomp.dylib <lib>` followed by `codesign -f -s -`).

### deal.II

```bash
git clone https://github.com/dealii/dealii.git
V=<env-dir>/view          # Spack view from above
OPENBLAS=/opt/homebrew/opt/openblas/lib/libopenblas.dylib

cmake -S dealii -B build-dealii \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HOME/dealii-install \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/mpicc \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/mpicxx \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DDEAL_II_ALLOW_AUTODETECTION=OFF \
  -DDEAL_II_COMPONENT_EXAMPLES=OFF \
  -DDEAL_II_COMPONENT_DOCUMENTATION=OFF \
  -DDEAL_II_WITH_MPI=ON \
  -DDEAL_II_WITH_LAPACK=ON -DLAPACK_LIBRARIES=$OPENBLAS -DBLAS_LIBRARIES=$OPENBLAS \
  -DDEAL_II_WITH_SCALAPACK=ON -DSCALAPACK_DIR=$V \
  -DDEAL_II_WITH_MUMPS=ON -DMUMPS_DIR=$V \
  -DDEAL_II_WITH_METIS=ON -DMETIS_DIR=$V \
  -DDEAL_II_WITH_P4EST=ON -DP4EST_DIR=$V \
  -DDEAL_II_WITH_TRILINOS=ON -DTRILINOS_DIR=$V \
  -DDEAL_II_WITH_KOKKOS=ON -DKOKKOS_DIR=$V \
  -DDEAL_II_WITH_BOOST=ON -DBOOST_DIR=$V \
  -DDEAL_II_WITH_TBB=ON -DTBB_DIR=$V \
  -DDEAL_II_WITH_ZLIB=ON -DZLIB_DIR=$V \
  -DDEAL_II_WITH_TASKFLOW=ON \
  -DDEAL_II_WITH_MUPARSER=ON \
  -DDEAL_II_WITH_UMFPACK=ON \
  -DOpenMP_ROOT=/opt/homebrew/opt/libomp \
  -DDEAL_II_WITH_PETSC=OFF -DDEAL_II_WITH_SLEPC=OFF

cmake --build build-dealii --target install -j 8
```

Check the configure summary: MPI, MUMPS, METIS, P4EST, TRILINOS (with Sacado),
BOOST, LAPACK and SCALAPACK must all be `ON`.

---

## 3. Building the code

Use an out-of-source build directory next to the repository:

```bash
mkdir build && cd build
cmake <path-to-repo> \
  -DDEAL_II_DIR=$HOME/dealii-install \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/libomp/lib -lomp"
make -j 4
```

Compiling the single large source file takes several minutes. The result is the
executable `build/nonlinear-poro-viscoelasticity`.

After editing the source, `make` in the build directory is enough.

---

## 4. Running

```bash
mpirun -np <ranks> <path>/nonlinear-poro-viscoelasticity <threads_per_rank> <parameter_file.prm>
```

- `<threads_per_rank>`: use **1** (pure MPI). That was the fastest setup on the
  Mac (and on the cluster).
- Paths in the `.prm` files (`input directory`, `output directory`) are
  **relative to the directory you start from**. The test files below assume
  you start from the repository root.
- The output directory **must exist** before the run; the code does not create it.

### Test case 1: generic compression relaxation

`parameter_files/parameters_quarter_cylinder_compression_relaxation.prm`

Quarter cylinder (R = 4 mm, H = 4.667 mm), glued top and bottom, drained
lateral surface. The top plate moves down 0.6 mm (~13 % compression) within
6 s and is held until 60 s.

```bash
mkdir -p out_compression_relaxation
mpirun -np 6 build/nonlinear-poro-viscoelasticity 1 parameter_files/parameters_quarter_cylinder_compression_relaxation.prm
```

Reference (Mac M-series, 6 ranks, MumpsDirect): 50 time steps, ~40 s. Peak
force −0.0300 N at 6 s, relaxing to −0.0128 N at 60 s (full cylinder, i.e.
`4e-6 * reac_cons [2]`, see below).

### Test case 2: experimental Am loading protocol

`parameter_files/parameters_quarter_cylinder_Am_experiment.prm`

Same geometry, loaded with the measured displacement history of the amygdala
(Am) rheometer test (mean of 12 specimens, ~960 s): cyclic compression/tension,
compression relaxation, tension relaxation, up to ±0.7 mm.

- `parameter_files/data/Am_load_path.csv`: load path read by the code
  (`time` in s, `displacement` in mm, one computed time step per row, 978 steps).
- `parameter_files/data/Am_experiment.csv`: measured data for comparison
  (`time`, `displacement`, `force` in N for the full cylinder; empty `force` =
  not measured).

```bash
mkdir -p out_Am_experiment
mpirun -np 6 build/nonlinear-poro-viscoelasticity 1 parameter_files/parameters_quarter_cylinder_Am_experiment.prm
```

Reference: 978 time steps, 8 min on 6 ranks with `BLRNewt` (measured);
expect roughly 2–3x longer with `MumpsDirect` (estimated). The simulated force follows the experiment closely
(relative L2 difference ~26 %, mainly at the start of the tension relaxation).

Both files use a visco-Ogden material (two-term Ogden equilibrium part plus one
Ogden Maxwell branch) with parameters from an inverse fit to the Am data.

---

## 5. Output

In the output directory:

| File | Content |
|---|---|
| `console-output.sol` | copy of the console log (mesh, Newton iterations, timings) |
| `data-for-gnuplot.sol` | one comma-separated row per time step; column names are in the last `#` header line |
| `solution.pvd`, `solution_*.pvtu`, `solution_*.vtu` | ParaView output every `Time step number output` steps (open `solution.pvd`) |

Most useful columns in `data-for-gnuplot.sol`:

- `Time`
- `P0[*]`, `P1[*]`: displacements (x, y, z) and pore pressure at the tracked
  points (P0: top centre (0, 0, H); P1: lateral surface at mid-height (R, 0, H/2),
  so `P1[0]` is the radial bulging)
- **`reac_cons [2]`: axial reaction force** (residual-based, the most accurate
  force; use this one, not `reaction [2]`).
  Units: for `brain_rheometer_cyclic_tension_compression_exp_quarter` (test case 2)
  it is already converted to **N for the full cylinder**. For the other geometry
  types (test case 1) it is in Pa·mm² for the quarter model; multiply by `4e-6`
  to get N for the full cylinder.
- `fluid flow`, `flow_cons`: fluid flux over the drained boundary
- `J_min_cell`: smallest volume ratio (check for excessive compression)

---

## 6. Parameter file: the most important settings

| Section / entry | Meaning |
|---|---|
| `geometry` → `Geometry type` | load case. `brain_rheometer_relaxation_tension_compression_quarter`: ramp + hold; `brain_rheometer_cyclic_tension_compression_exp_quarter`: displacement history from a CSV file |
| `geometry` → `Global refinement` | mesh level (2 = default, 3 = finer; each level ~8x more DoFs) |
| `geometry` → `Edge refinement` | extra refinement at the glued edges |
| `geometry` → `height`, `radius` | specimen size in mm |
| `testing_device` → `Load value` | final displacement in mm (ramp + hold case; negative = compression) |
| `testing_device` → `input files` | load-path CSV (with `global` → `input directory`) |
| `Material properties` → `material` | `Neo-Hooke`, `Ogden`, `visco-Ogden`, `visco2-Ogden` |
| `Material properties` → `lambda`, `mu*`, `alpha*`, `viscosity_*`, `initial intrinsic permeability` | material parameters (Pa, Pa·s, mm²) |
| `Nonlinear solver` → `Linear solver` | `MumpsDirect` (robust, default in the test files) or `BLRNewt` (faster: low-rank MUMPS as preconditioner for FGMRES) |
| `Time` → `End time`, `End load time`, `Time step size` | time stepping for the ramp + hold case (ignored if a load-path CSV gives the time steps) |
| `Output parameters` → `Output files` | `none`, `solution`, `bcs`, `all` |
