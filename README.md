# stl2step

Convert STL meshes into STEP (AP214) B-rep files that import into SolidWorks,
Fusion, FreeCAD, Onshape and other CAD tools.

## What it does
- Reads binary and ASCII STL
- Welds vertices, removes degenerate/duplicate triangles, fixes inverted winding
- Merges coplanar neighbouring triangles into real planar faces (with holes)
- Writes a solid (watertight mesh) or a surface body (open mesh)

## Limitation
An STL contains only triangles. Curved surfaces (holes, fillets, cylinders)
remain faceted in the STEP. You get a usable solid, not a parametric feature
tree. See the roadmap below.

## Build
Requires [Rust](https://rustup.rs) (on Windows, with the MSVC build tools).

```
cargo build --release                          # GUI exe
cargo build --release --no-default-features    # console-only exe
cargo test --no-default-features               # unit tests
```
Output: `target/release/stl2step(.exe)`.

## Use
GUI: run the exe, drag an `.stl` onto the window (or the exe icon), Convert.

Console build:
```
stl2step part.stl -o part.step [--scale 25.4] [--no-merge] [--angle-tol 0.1]
```
`--scale 25.4` for inch-based STLs (output is always millimetres).

## Releases
Pushing a tag like `v1.0.1` makes GitHub Actions build `stl2step.exe` and
attach it to a GitHub Release automatically.

## Layout
- `src/converter.rs` conversion core (std only)
- `src/gui.rs`, `src/cli.rs` front ends
- `python/stl2step.py` original reference implementation (stdlib only)

## Roadmap
- [ ] Cylinder detection (true round holes/bosses)
- [ ] Merge collinear edges on merged faces
- [ ] Mesh hole-filling for non-watertight STLs
- [ ] Sphere / cone / fillet fitting
- [ ] Multi-body STL (split disconnected shells)

## License
MIT. See `LICENSE`.
