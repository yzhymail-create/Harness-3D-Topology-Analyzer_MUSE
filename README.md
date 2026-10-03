# Harness 3D Topology Analyzer (MUSE working version)

One-click: CATIA-exported STEP → 6-sheet wiring-harness FROMTO Excel report.

Validated end-to-end on X1X-DRD (22 solids: 14 curves, 7 connectors, 21/21 paths OK,
zero `#S` residue, all length-conservation checks pass).
IP-Harness (189 solids) name mapping validated; full run pending.

## Files

| File | Role |
|---|---|
| `harness_3d/fromto_report.py` | Main entry: STEP → 6-sheet XLSX (`Meta/Connectors/Curves/FromToPaths/FromToLatest/Diagnostics`) |
| `harness_3d/harness_topology.py` | Contact-based topology extraction (tubes, connectors, clamps, stations, segments, nodes, runs) |
| `harness_3d/map_step_names.py` | STEP text → real instance names (`S{i}` mapping, standalone usable) |
| `harness_3d/tests/test_fromto.py` | 23 FROMTO report tests |
| `harness_3d/tests/test_topology_logic.py` | 50 topology logic tests |
| `harness_3d/requirements.txt` | Python 3.8 compatible deps |
| `harness_3d/requirements-py312.txt` | Python 3.12 deps |

## Usage

```bash
pip install -r harness_3d/requirements.txt   # or requirements-py312.txt
python harness_3d/fromto_report.py <input.stp> -o <output.xlsx>
# optional: --original <old_fromto.xlsx>  (adds OrigLength comparison column)
```

Tests:

```bash
python harness_3d/tests/test_topology_logic.py
python harness_3d/tests/test_fromto.py
```

## Key design decisions (validated, do not "simplify" away)

### 1. Name mapping — `S{i}` = representation item order, NOT definition order

- OCCT transfer order follows the `ADVANCED_BREP_SHAPE_REPRESENTATION` item list,
  **not** the order `MANIFOLD_SOLID_BREP` entities appear in the STEP text.
- Both `MANIFOLD_SOLID_BREP` **and** `BREP_WITH_VOIDS` must be parsed
  (IP-Harness has 3 `BREP_WITH_VOIDS`; missing them caused a 186/189 mismatch).
- `\X2\...\X0\` escape sequences are decoded to real CJK characters in-code.
- The tool self-checks: parsed name count must equal OCP solid count, otherwise it
  aborts **before** the hours-long topology pass.
- Do NOT match names via OCP `TransientProcess.Find(entity)` on sub-entities
  (returns empty) or via `model.Value(n)` numbering (≠ STEP `#id`).
  Do NOT match via vertex-coordinate lookup — it silently fails on `BREP_WITH_VOIDS`.

### 2. Positioning rule — derive from the TUBE, never from the device itself

A clamp/connector's reported coordinate is **not** its volume center, bbox center,
or any feature point on its own geometry. It is derived from the tube it contacts:

- Device touches a tube **end** (within 5 mm of the end) → coordinate = that
  tube segment's **end-face center**.
- Clamp touches a tube **middle** → coordinate = the **nearest** tube segment's
  end-face center (topology internally splits at the spine projection point;
  the report caliber is uniformly end-face centers).

The device's own solid is used only for contact detection ("does it touch?"),
never for positioning. Volume-center-based approaches are wrong in principle.

### 3. Thresholds (mm, validated against CATIA VBA ground truth)

| Constant | Value | Meaning |
|---|---|---|
| `CONTACT_TOL` | 2.0 | solid contact decision distance |
| `END_TOL` | 5.0 | support point within this of tube end ⇒ end contact |
| `STATION_END_TOL` | 5.0 | tap projection within this of main-tube end ⇒ treated as end joint |
| station merge | 10.0 | arc-length + spatial dual condition |

### 4. Contact performance — 3-layer filter, `BRepExtrema` only at the last layer

`min_contact()` in `harness_topology.py`:

1. Bounding-sphere prefilter (cheap but loose — nearly useless for long tubes).
2. **Exact** bounding-box distance via `Bnd_Box.Distance` (tight and provably safe:
   box distance > `CONTACT_TOL` ⇒ true distance > `CONTACT_TOL`, so skipping
   cannot change the result).
3. `BRepExtrema_DistShapeShape` only for pairs passing (1) and (2).

Layer 2 was verified to produce **byte-identical** report output on X1X
(all 5 data sheets, cell-by-cell) while eliminating the bulk of exact distance calls.

### 5. Report format — 6 sheets, fixed

`Meta / Connectors / Curves / FromToPaths / FromToLatest / Diagnostics`.
`FromToPaths.Status` ∈ {`OK`, `NOT_CONNECTED`}; broken pairs carry empty
Length/PathMarker with gap details in `Diagnostics.BrokenDetails`.
`PathMarker` = JSON `[connector, curve, (clamp, curve)*, connector]`.
Path length = sum of curve lengths (conservation checked per row).

### 6. Stable coding

`SEG / CON / TIE / BN / N`. Clamp code is `TIE` (not `CLP`).
