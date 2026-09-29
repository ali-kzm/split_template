# Offline Cut-Cell Template Generator

Standalone C++20 generator for validated, quad-dominant mixed triangle/quadrilateral cut-cell templates on the closed reference square `[-1,1] x [-1,1]`. It is intentionally independent of RheoFEM.

## Supported canonical families

Corner order is counterclockwise: `0=(-1,-1), 1=(1,-1), 2=(1,1), 3=(-1,1)`.

- case 2: `-+++`, interface crosses edges `(0,1)` and `(3,0)`
- case 3: `--++`, interface crosses edges `(1,2)` and `(3,0)`
- case 6: `0-++`, interface connects corner 0 to edge `(1,2)`
- case 11: `0+0-`, fixed diagonal from corner 0 to corner 2

Cases 5 and 8 are intentionally outside this generator. Negative `phi` is the inside phase. `phi` is the signed distance to the supporting line. `phi_zero_tol` is an absolute tolerance in reference-coordinate signed-distance units and is independent of geometric/quality tolerances.

## Dependencies

Tested design targets:

- C++20 compiler (GCC 12+ or Clang 15+)
- CMake 3.20+
- CGAL 5.5+ (GPL/LGPL components; this project uses the triangulation packages)
- Boost.Graph 1.74+ (Boost Software License 1.0)
- libjpeg-turbo 2.1+ (BSD-style/IJG/zlib licenses as distributed upstream)

Ubuntu/Debian:

```sh
sudo apt update
sudo apt install -y build-essential cmake ninja-build libcgal-dev libboost-graph-dev libjpeg-dev
```

On Ubuntu, `libjpeg-dev` resolves to the system JPEG development implementation, normally libjpeg-turbo. CMake must find a real JPEG library; the program writes genuine JPEG bitstreams through the libjpeg API.

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
  --min-triangle-quality 0.05 \
  --min-quad-quality 0.2 \
  --max-nodes 512 \
  --max-attempts-per-template 100 \
  --max-smoothing-passes 50 \
  --image-size 1024
```

The generator uses one constrained triangulation for both phases, maximum-weight matching for admissible triangle pairs, constrained Laplacian smoothing with backtracking, validation, enrichment, phase-specific pressure numbering, then atomic publication of each `.dat`/JPEG pair.

Quad dominance is an objective for the current admissible pairing graph, not a claim of globally optimal quadrilateral meshing. Valid unmatched triangles are retained. A finite template set does not cover arbitrarily degenerate cuts.

See [docs/DATASET_FORMAT.md](docs/DATASET_FORMAT.md) for the versioned data contract.

## Scope boundary

This repository does **not** implement runtime template matching, runtime coordinate adjustment, post-adjustment validation in RheoFEM, global FEM numbering, MLS coefficients, or RheoFEM integration.
