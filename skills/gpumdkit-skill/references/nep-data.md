# NEP Training and Prediction Data

Use this reference to create and audit `train.xyz` and `test.xyz`.

## Contents

- Reference scope
- Configuration format
- Targets by model type
- Units and conventions
- Dataset audit and split

## Reference scope

This file contains the portable `train.xyz` and `test.xyz` specification. Both files use multi-frame extended XYZ; each configuration occupies `N + 2` lines.

## Configuration format

Line 1:

```text
<number_of_atoms>
```

Line 2 contains header fields. `lattice`, `energy`, and `properties` are parser-required for every configuration, even when a non-potential model later ignores energy and force-related targets:

| Field | Meaning |
|---|---|
| `lattice="ax ay az bx by bz cx cy cz"` | Mandatory cell vectors, Angstrom |
| `energy=<total_energy>` | Mandatory total-energy field, eV; ignored for dipole/polarizability training |
| `virial="9 components"` | Optional total-cell virial tensor, eV |
| `stress="9 components"` | Optional stress tensor, eV/Angstrom^3; ignored if virial is present |
| `weight=<relative_weight>` | Optional configuration weight |
| `properties=...` | Mandatory per-atom column schema |

Relevant per-atom properties:

```text
species:S:1
pos:R:3
force:R:3        # or forces:R:3
bec:R:9          # optional Born effective charge target
```

Example schema only:

```text
properties=species:S:1:pos:R:3:force:R:3
```

Do not reuse example energies, cells, or structures. Confirm each atom row matches the declared schema and line-1 atom count.

## Targets by model type

| Model | Required/used targets |
|---|---|
| Potential (`model_type 0`) | Energy and forces; virial/stress when included in the intended loss |
| Dipole (`model_type 1`) | Header `dipole="dx dy dz"`; energy/virial/stress/force ignored |
| Polarizability (`model_type 2`) | Header `pol="9 components"`; energy/virial/stress/force/dipole ignored |
| qNEP potential | Potential targets; target charges are not required; optional `bec:R:9` for BEC fitting |

For dipole and polarizability, confirm whether targets are global or per-atom and set `atomic_v` consistently. Do not mix incompatible target conventions.

## Units and conventions

| Quantity | Unit/convention |
|---|---|
| Length/position | Angstrom |
| Energy | eV, total per cell in input |
| Force | eV/Angstrom |
| Virial | eV, total per cell in input |
| Stress | eV/Angstrom^3 in input |
| BEC | elementary charge |
| Dipole/polarizability | User-selected consistent unit; record it |

Virial input uses positive values for compression and negative for tension. Stress input uses the opposite sign convention: positive for tension and negative for compression. If both are present, virial is used.

NEP training assumes periodic boundaries in all directions. For slabs, clusters, or vacuum systems, ask the user how the periodic cell and reference calculation are intended to represent the system.

The code may internally replicate cells thinner than twice the radial cutoff. Verify that this is physically appropriate.

Training uses single precision. Reference energies below `-100 eV/atom` can lose accuracy. Report this limitation and ask how reference energies should be handled; do not silently shift them.

## Dataset audit and split

- Verify species symbols and exact `type` order.
- Verify energy is total per cell, not per atom; verify virial/stress sign and units.
- Confirm force columns and atom ordering for every frame.
- Detect NaN/Inf, malformed cells, extreme distances, duplicated configurations, and inconsistent labels.
- Inspect energy/force/virial ranges by composition and configuration class.
- Keep trajectory-neighbor configurations from leaking across train/test when independence matters.
- Preserve raw data and record filters, weights, exclusions, and random seed.
- Ensure every intended species, phase, composition, strain, defect, surface, temperature, and short-range regime is represented according to the user's model scope.
- Do not invent a train/test fraction or minimum dataset size.

Use `references/format-conversion.md`, `analyzers.md`, and `sampling.md` for GPUMDkit commands after confirming the transformation criteria.
