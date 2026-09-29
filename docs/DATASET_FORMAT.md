# PVMLS template dataset format

Schema version: **2.0.0**

Each accepted template is published as `case_<id>/temp_<NN>.dat` plus a matching genuine JPEG preview. `manifest.json` records accepted templates and failed generation requests.

The text file is line-oriented and human readable. Floating-point values are serialized with 17 significant decimal digits.

## Sections

- `PVMLS_TEMPLATE_DATA 2.0.0`
- `META`: case, template index, seed and complete generation settings
- `CORNERS`: canonical corner ID, coordinates, raw signed-distance phi and classified sign
- `INTERFACE`: endpoints, explicit prescribed endpoint `phi=0`, supporting-line data and edge sampling parameters
- `NODES`: phase-owned geometric nodes, including enriched nodes
- `PRESSURE_RECORDS`: local pressure ID -> geometric node ID + phase
- `TRIANGLES`: phase, quality, P1 vertex connectivity, P2 connectivity
- `QUADS`: phase, quality, Q1 vertex connectivity, Q2 connectivity
- `PRESSURE_CONNECTIVITY`: pressure P1/Q1 connectivity per cell
- `BOUNDARY_EDGES`: ordered primary-mesh square-boundary segments
- `INTERFACE_EDGES_NEGATIVE`: ordered negative-side primary interface segments
- `INTERFACE_EDGES_POSITIVE`: ordered positive-side primary interface segments
- `METRICS`: min/mean triangle and quad quality, quad count fraction, quad area fraction, per-phase areas, total primary element count, negative/positive element counts, target element count, and hard maximum
- `VALIDATION OK`
- `END`

## Node fields

`id ksi eta phi sign owner_phase constraint parent_edge parent_segment parameter primary`

Constraint is one of `corner`, `boundary`, `interface`, or `interior`. `owner_phase` is `-1` or `+1` and identifies which independently meshed phase owns the geometric node. `parent_edge=-1` means no square edge. For interface nodes, `parent_segment=0` denotes the negative-side interface chain and `parent_segment=1` denotes the positive-side chain.

The two phases use the same straight geometric interface and the same prescribed endpoints, but their interior interface discretizations are deliberately **nonmatching**. Internal interface geometric node IDs are not shared across phases; one phase may contain more interface nodes than the other.

Primary nodes are optimized first. Enriched nodes are then created as exact edge midpoints and quadrilateral bilinear centers. Shared edge midpoints are globally deduplicated within one template.

## Connectivity conventions

- Triangle P1: `v0 v1 v2`, counterclockwise.
- Triangle P2: `v0 v1 v2 m01 m12 m20`.
- Quadrilateral Q1: `v0 v1 v2 v3`, counterclockwise.
- Quadrilateral Q2: `v0 v1 v2 v3 m01 m12 m23 m30 center`.

Pressure connectivity uses P1 for triangles and Q1 for quadrilaterals. Because geometric nodes are phase-owned in schema 2.0, nonmatching interface nodes naturally have separate geometric and pressure identities on each side.

## Quality definitions

Triangle quality:

`4*sqrt(3)*area / (sum of squared edge lengths)`.

Quadrilateral quality is the minimum signed corner scaled Jacobian for counterclockwise vertices. At corner `i`, the two vectors point from the corner to its next and previous vertices, and the signed scaled Jacobian is their 2-D cross product divided by the product of their lengths.

## Phi classification

With positive finite `phi_zero_tol`:

- `phi < -phi_zero_tol`: negative
- `phi > +phi_zero_tol`: positive
- otherwise: zero

Raw signed-distance values are retained for ordinary nodes. Prescribed interface nodes store exactly zero.
