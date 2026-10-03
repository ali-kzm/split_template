# Offline Cut-Cell Template Generator

Standalone C++20 generator for validated, quad-dominant mixed triangle/quadrilateral cut-cell templates on the closed reference square `[-1,1] x [-1,1]`. It is intentionally independent of RheoFEM.

## Supported canonical families

Corner order is counterclockwise: `0=(-1,-1), 1=(1,-1), 2=(1,1), 3=(-1,1)`.

- case 2: `-+++`, interface crosses edges `(0,1)` and `(3,0)`
- case 3: `--++`, interface crosses edges `(1,2)` and `(3,0)`
- case 6: `0-++`, interface connects corner 0 to edge `(1,2)`
- case 11: `0+0-`, fixed diagonal from corner 0 to corner 2

Cases 5 and 8 are intentionally outside this generator. Negative `phi` is the inside phase. `phi` is the signed distance to the supporting line. `phi_zero_tol` is an absolute tolerance in reference-coordinate signed-distance units and is independent of geometric/quality tolerances.

## Dependencies and installation

### Libraries used by this project

The project directly uses:

| Dependency | Purpose | CMake target / lookup |
| --- | --- | --- |
| C++20 compiler | Standard language/runtime support | C++20 |
| CMake 3.20+ | Configure and generate the build | — |
| CGAL 5.5+ | Robust predicates and constrained Delaunay triangulation | `CGAL::CGAL` |
| Boost.Graph 1.74+ | Maximum-weight matching for triangle pairing | `Boost::graph` |
| libjpeg-turbo / libjpeg API | Genuine JPEG preview encoding and decoding in tests | `JPEG::JPEG` |

CGAL also uses **GMP** and **MPFR** transitively. Package managers install these automatically with CGAL, so they normally do not need to be installed separately.

This project does **not** require Qt, Eigen, OpenGL, RheoFEM, MPI, TBB, or OpenCL.

Relevant licenses:

- CGAL: package components are distributed under GPL/LGPL and other compatible component licenses; check the CGAL package/module used by your redistribution.
- Boost.Graph: Boost Software License 1.0.
- libjpeg-turbo: BSD-3-Clause and IJG licensing as distributed upstream.
- GMP: LGPL/GPL dual licensing.
- MPFR: LGPL.

### Ubuntu / Debian

Ubuntu 22.04/24.04 and recent Debian releases can install everything from the system package manager:

```sh
sudo apt update
sudo apt install -y \
  build-essential \
  cmake \
  ninja-build \
  libcgal-dev \
  libboost-graph-dev \
  libjpeg-dev
```

`libcgal-dev` installs the required CGAL development files and pulls in GMP/MPFR dependencies. On Ubuntu, `libjpeg-dev` is provided by the system JPEG development implementation, normally libjpeg-turbo.

Configure, build, and test:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

### Windows

The recommended Windows setup is **Visual Studio 2022 + vcpkg**.

1. Install **Visual Studio 2022** or **Visual Studio 2022 Build Tools** with the **Desktop development with C++** workload. Make sure MSVC, the Windows SDK, and CMake support are enabled.

2. Install Git if it is not already available, then install vcpkg from PowerShell:

```powershell
git clone https://github.com/microsoft/vcpkg C:\src\vcpkg
C:\src\vcpkg\bootstrap-vcpkg.bat
```

3. Install the project libraries:

```powershell
C:\src\vcpkg\vcpkg.exe install cgal:x64-windows boost-graph:x64-windows libjpeg-turbo:x64-windows
```

vcpkg installs CGAL's required GMP/MPFR and Boost dependencies automatically.

4. Configure the project using the vcpkg CMake toolchain:

```powershell
cmake -S . -B build -A x64 ^
  -DCMAKE_TOOLCHAIN_FILE=C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake
```

If running the command directly in PowerShell instead of `cmd.exe`, use one line:

```powershell
cmake -S . -B build -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake
```

5. Build and run the tests:

```powershell
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The executable is normally produced under:

```text
build/Release/pvmls-template-generator.exe
```

On Linux/Ninja builds it is normally:

```text
build/pvmls-template-generator
```

## Build and test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Generate the default dataset

```sh
./build/pvmls-template-generator
```

This requests 15 accepted templates for each of cases 2, 3, 6, and 11 in `./pvmls_templates`. Existing generated content is not overwritten unless `--overwrite` is passed.

Example:

```sh
./build/pvmls-template-generator \
  --output ./pvmls_templates \
  --templates-per-case 15 \
  --cases 2,3,6,11 \
  --seed 24301 \
  --phi-zero-tol 1e-10 \
  --min-edge-fraction 1e-3 \
  --target-edge-length 0.3 \
  --target-elements 10 \
  --max-elements 20 \
  --triangle-aspect-threshold 3 \
  --quad-aspect-threshold 3 \
  --pentagon-small-edge-fraction 0.15 \
  --min-triangle-quality 0.05 \
  --min-quad-quality 0.2 \
  --max-nodes 512 \
  --max-attempts-per-template 100 \
  --max-smoothing-passes 50 \
  --image-size 1024
```

The generator meshes the negative and positive clipped phase polygons independently, then combines them into one template. Both phase boundaries lie on the same straight interface and share the prescribed interface endpoints, but **interior interface nodes are intentionally nonconforming**: one side may have a different number of interface nodes from the other.

### Polygon-specific meshing strategy

The clipped phase polygon is classified before meshing. Only convex triangles, quadrilaterals and pentagons occur for the supported straight cuts.

- **Triangle:** aspect ratio is `longest_edge^2 / (2*area)`. If it is below `--triangle-aspect-threshold` (default **3**), the polygon is kept as one triangle. A thin triangle is peeled from its wide end into quadrilateral strips; one geometrically similar narrow triangle is deliberately retained at the tip.
- **Quadrilateral:** regular quads (`aspect <= --quad-aspect-threshold`, default **3**) are split into **2 or 4 quads** according to size. Thin/distorted quads use a structured bilinear grid with more subdivisions along the long direction, bounded by the phase and global element budgets.
- **Pentagon:** all five possible triangle+quad ear diagonals are assessed. If a boundary edge is very small (below `--pentagon-small-edge-fraction`, default **0.15** of the longest edge), the generator first tries to isolate that short edge in a local ear triangle while keeping the remaining quadrilateral valid. Clean non-tiny pentagons may also use the best valid 2-cell split. If no safe ear split exists, the fallback is a **5-quad center/mid-edge decomposition**.

The polygon strategy directly constructs the final primary triangles/quads; triangle pairing is no longer part of the active generation path. Constrained smoothing, validation, enrichment, pressure numbering and atomic publication still run afterward.

Mesh size remains bounded by `--max-elements` (default **20**). `--target-elements` is now a sizing hint mainly used for distorted quadrilateral subdivision; regular triangles intentionally ignore it and remain single cells. `--target-edge-length` remains a legacy compatibility option.

**Case 11 note:** the fixed diagonal creates two regular triangles, both below the default aspect threshold. Therefore the new rules produce one canonical two-triangle template. The generator emits one unique case-11 template even when `--templates-per-case` is larger.

See [docs/DATASET_FORMAT.md](docs/DATASET_FORMAT.md) for the versioned data contract.

## Scope boundary

This repository does **not** implement runtime template matching, runtime coordinate adjustment, post-adjustment validation in RheoFEM, global FEM numbering, MLS coefficients, or RheoFEM integration.
