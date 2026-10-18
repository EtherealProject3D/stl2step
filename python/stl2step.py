#!/usr/bin/env python3
"""
stl2step - convert STL meshes into STEP (AP214) B-rep solids.

* Pure Python 3 standard library (tkinter is used for the GUI if available).
* Reads binary and ASCII STL.
* Welds duplicate vertices, drops degenerate triangles, fixes inverted winding.
* Merges coplanar neighbouring triangles into single planar faces (with holes),
  so flat areas arrive in SolidWorks as real faces you can select, sketch on,
  extrude from, etc. instead of thousands of little triangles.
* Writes a MANIFOLD_SOLID_BREP (watertight meshes) or a shell-based surface
  model (open meshes).

Usage
    python stl2step.py                       -> opens the GUI
    python stl2step.py part.stl              -> writes part.step
    python stl2step.py part.stl -o out.step --scale 25.4 --no-merge

IMPORTANT: an STL has no curves, only triangles. Curved surfaces (holes,
fillets, cylinders) will stay faceted in the STEP file. Converting those back
into true cylinders/fillets/features is "feature recognition" / reverse
engineering and needs dedicated CAD tools (see README notes in the GUI).
"""
import argparse
import datetime
import math
import os
import re
import struct
import sys
from collections import defaultdict

# --------------------------------------------------------------------------
# small vector helpers
# --------------------------------------------------------------------------

def sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])
def norm(a): return math.sqrt(dot(a, a))
def unit(a):
    n = norm(a)
    return (a[0] / n, a[1] / n, a[2] / n) if n > 0 else (0.0, 0.0, 1.0)


# --------------------------------------------------------------------------
# STL reading
# --------------------------------------------------------------------------

_VERTEX_RE = re.compile(
    rb"vertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)")


def read_stl(path):
    """Return list of triangles ((x,y,z),(x,y,z),(x,y,z))."""
    with open(path, "rb") as f:
        data = f.read()
    if len(data) < 15:
        raise ValueError("File is too small to be an STL.")

    tris = []
    is_binary = False
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if 84 + n * 50 == len(data):
            is_binary = True
        elif not data[:5].lower() == b"solid":
            # binary with a wrong/truncated count
            is_binary = True
            n = (len(data) - 84) // 50
    if is_binary:
        n = min(struct.unpack_from("<I", data, 80)[0], (len(data) - 84) // 50)
        for rec in struct.iter_unpack("<12fH", data[84:84 + n * 50]):
            tris.append(((rec[3], rec[4], rec[5]),
                         (rec[6], rec[7], rec[8]),
                         (rec[9], rec[10], rec[11])))
    else:
        vs = [(float(a), float(b), float(c))
              for a, b, c in _VERTEX_RE.findall(data)]
        if len(vs) % 3:
            raise ValueError("ASCII STL has a vertex count not divisible by 3.")
        for i in range(0, len(vs), 3):
            tris.append((vs[i], vs[i + 1], vs[i + 2]))
    if not tris:
        raise ValueError("No triangles found in file.")
    return tris


# --------------------------------------------------------------------------
# mesh cleanup
# --------------------------------------------------------------------------

class Mesh:
    def __init__(self):
        self.verts = []   # list of (x,y,z)
        self.faces = []   # list of (a,b,c) indices, CCW seen from outside
        self.closed = False
        self.warnings = []


def build_mesh(tris, weld_rel=1e-5, scale=1.0, log=print):
    m = Mesh()
    if scale != 1.0:
        tris = [tuple((p[0] * scale, p[1] * scale, p[2] * scale) for p in t)
                for t in tris]

    xs = [p[0] for t in tris for p in t]
    ys = [p[1] for t in tris for p in t]
    zs = [p[2] for t in tris for p in t]
    diag = math.sqrt((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2 +
                     (max(zs) - min(zs)) ** 2) or 1.0
    m.diag = diag
    tol = diag * weld_rel
    inv = 1.0 / tol

    vid = {}
    seen = set()
    degenerate = dup = 0
    min_area2 = (1e-12 * diag * diag) ** 2
    for t in tris:
        idx = []
        for p in t:
            key = (round(p[0] * inv), round(p[1] * inv), round(p[2] * inv))
            i = vid.get(key)
            if i is None:
                i = len(m.verts)
                vid[key] = i
                m.verts.append(p)
            idx.append(i)
        a, b, c = idx
        if a == b or b == c or a == c:
            degenerate += 1
            continue
        cr = cross(sub(m.verts[b], m.verts[a]), sub(m.verts[c], m.verts[a]))
        if dot(cr, cr) < min_area2 * 4:
            degenerate += 1
            continue
        key = (a, b, c)
        rot = min((a, b, c), (b, c, a), (c, a, b))
        if rot in seen:
            dup += 1
            continue
        seen.add(rot)
        m.faces.append((a, b, c))

    if degenerate:
        m.warnings.append(f"Removed {degenerate} degenerate triangle(s).")
    if dup:
        m.warnings.append(f"Removed {dup} duplicate triangle(s).")

    # topology check
    dcount = defaultdict(int)
    for a, b, c in m.faces:
        dcount[(a, b)] += 1
        dcount[(b, c)] += 1
        dcount[(c, a)] += 1
    m.nonmanifold = any(v > 1 for v in dcount.values())
    m.closed = (not m.nonmanifold) and all(
        dcount.get((b, a), 0) == 1 for (a, b) in dcount)

    if m.closed:
        vol = 0.0
        for a, b, c in m.faces:
            vol += dot(m.verts[a], cross(m.verts[b], m.verts[c]))
        if vol < 0:
            m.faces = [(a, c, b) for a, b, c in m.faces]
            m.warnings.append("Triangle winding was inverted; flipped it.")
    else:
        if m.nonmanifold:
            m.warnings.append(
                "Mesh has non-manifold edges. Output will be a surface "
                "model, not a solid.")
        else:
            m.warnings.append(
                "Mesh is not watertight (open edges). Output will be a "
                "surface model, not a solid. Repair holes in the STL if you "
                "need a solid.")
    return m


# --------------------------------------------------------------------------
# planar face merging
# --------------------------------------------------------------------------

def _tri_normal(verts, f):
    a, b, c = (verts[i] for i in f)
    return unit(cross(sub(b, a), sub(c, a)))


def merge_planar(m, ang_tol_deg=0.1, dist_rel=1e-4):
    """
    Group coplanar, edge-connected triangles.
    Returns list of (normal, outer_loop, [hole_loops]) where loops are lists
    of vertex indices (outer CCW, holes CW about the normal).
    Regions that cannot be turned into clean loops fall back to single
    triangles.
    """
    verts, faces = m.verts, m.faces
    n = len(faces)
    normals = [_tri_normal(verts, f) for f in faces]
    cos_tol = math.cos(math.radians(ang_tol_deg))
    dist_tol = dist_rel * m.diag

    emap = {}
    for t, (a, b, c) in enumerate(faces):
        emap[(a, b)] = t
        emap[(b, c)] = t
        emap[(c, a)] = t

    visited = [False] * n
    out = []
    merged_tris = 0
    for seed in range(n):
        if visited[seed]:
            continue
        n0 = normals[seed]
        d0 = dot(n0, verts[faces[seed][0]])
        visited[seed] = True
        region = [seed]
        stack = [seed]
        while stack:
            t = stack.pop()
            a, b, c = faces[t]
            for u, v in ((a, b), (b, c), (c, a)):
                nb = emap.get((v, u))
                if nb is None or visited[nb]:
                    continue
                if dot(normals[nb], n0) < cos_tol:
                    continue
                if any(abs(dot(n0, verts[i]) - d0) > dist_tol
                       for i in faces[nb]):
                    continue
                visited[nb] = True
                region.append(nb)
                stack.append(nb)

        loops = None
        if len(region) > 1:
            loops = _region_loops(region, faces, verts, n0)
        if loops is None:
            for t in region:
                a, b, c = faces[t]
                out.append((normals[t], [a, b, c], []))
        else:
            merged_tris += len(region)
            out.append((n0, loops[0], loops[1:]))
    return out, merged_tris


def _region_loops(region, faces, verts, n0):
    dset = set()
    for t in region:
        a, b, c = faces[t]
        dset.add((a, b)); dset.add((b, c)); dset.add((c, a))
    nxt = defaultdict(list)
    for (a, b) in dset:
        if (b, a) not in dset:
            nxt[a].append(b)
    if any(len(v) != 1 for v in nxt.values()):
        return None  # pinch point; fall back
    loops, seen = [], set()
    for s in nxt:
        if s in seen:
            continue
        loop = [s]
        seen.add(s)
        cur = nxt[s][0]
        while cur != s:
            if cur in seen or cur not in nxt:
                return None
            loop.append(cur)
            seen.add(cur)
            cur = nxt[cur][0]
        loops.append(loop)

    def area(loop):
        tot = 0.0
        for i in range(len(loop)):
            p, q = verts[loop[i]], verts[loop[(i + 1) % len(loop)]]
            tot += dot(cross(p, q), n0)
        return 0.5 * tot

    areas = [area(l) for l in loops]
    pos = [i for i, a in enumerate(areas) if a > 0]
    if len(pos) != 1:
        return None
    o = pos[0]
    return [loops[o]] + [l for i, l in enumerate(loops) if i != o]


# --------------------------------------------------------------------------
# STEP writer
# --------------------------------------------------------------------------

def _r(x):
    s = repr(float(x))
    if "e" in s or "E" in s:
        mant, exp = s.lower().split("e")
        if "." not in mant:
            mant += "."
        return f"{mant}E{int(exp):+d}".replace("+", "")
    if "." not in s:
        s += "."
    return s


def _q(s):
    return s.replace("'", "''")


class _Step:
    def __init__(self):
        self.lines = []

    def add(self, body):
        self.lines.append(body)
        return len(self.lines)  # id = index (1-based)

    def ref(self, i): return f"#{i}"


def write_step(m, faces_geo, path, name="part"):
    S = _Step()
    add = S.add
    verts = m.verts
    name = _q(name)

    # --- units / context
    lu = add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))")
    au = add("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))")
    su = add("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())")
    unc = add(f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#{lu},"
              "'distance_accuracy_value','confusion accuracy')")
    ctx = add("(GEOMETRIC_REPRESENTATION_CONTEXT(3) "
              f"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{unc})) "
              f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{lu},#{au},#{su})) "
              "REPRESENTATION_CONTEXT('Context3D','3D Context with UNIT "
              "and UNCERTAINTY'))")

    # --- vertices
    pt, vx = [], []
    for p in verts:
        i = add(f"CARTESIAN_POINT('',({_r(p[0])},{_r(p[1])},{_r(p[2])}))")
        pt.append(i)
        vx.append(add(f"VERTEX_POINT('',#{i})"))

    # --- edges (one EDGE_CURVE per undirected edge)
    edges = {}

    def edge_curve(a, b):
        key = (a, b) if a < b else (b, a)
        ec = edges.get(key)
        if ec is None:
            lo, hi = key
            d = sub(verts[hi], verts[lo])
            L = norm(d)
            u = unit(d)
            di = add(f"DIRECTION('',({_r(u[0])},{_r(u[1])},{_r(u[2])}))")
            ve = add(f"VECTOR('',#{di},{_r(L)})")
            ln = add(f"LINE('',#{pt[lo]},#{ve})")
            ec = add(f"EDGE_CURVE('',#{vx[lo]},#{vx[hi]},#{ln},.T.)")
            edges[key] = ec
        return ec

    def loop_id(loop):
        oes = []
        k = len(loop)
        for i in range(k):
            a, b = loop[i], loop[(i + 1) % k]
            ec = edge_curve(a, b)
            oes.append(add(f"ORIENTED_EDGE('',*,*,#{ec},"
                           f"{'.T.' if a < b else '.F.'})"))
        return add("EDGE_LOOP('',(" + ",".join(f"#{o}" for o in oes) + "))")

    # --- faces
    face_ids = []
    for normal, outer, holes in faces_geo:
        o = add(f"CARTESIAN_POINT('',({_r(verts[outer[0]][0])},"
                f"{_r(verts[outer[0]][1])},{_r(verts[outer[0]][2])}))")
        # reference direction: in-plane axis
        axis = min(range(3), key=lambda i: abs(normal[i]))
        helper = [0.0, 0.0, 0.0]
        helper[axis] = 1.0
        ref = unit(cross(normal, helper))
        zd = add(f"DIRECTION('',({_r(normal[0])},{_r(normal[1])},"
                 f"{_r(normal[2])}))")
        xd = add(f"DIRECTION('',({_r(ref[0])},{_r(ref[1])},{_r(ref[2])}))")
        ax = add(f"AXIS2_PLACEMENT_3D('',#{o},#{zd},#{xd})")
        pl = add(f"PLANE('',#{ax})")
        bounds = [add(f"FACE_OUTER_BOUND('',#{loop_id(outer)},.T.)")]
        for h in holes:
            bounds.append(add(f"FACE_BOUND('',#{loop_id(h)},.T.)"))
        face_ids.append(add("ADVANCED_FACE('',("
                           + ",".join(f"#{b}" for b in bounds)
                           + f"),#{pl},.T.)"))

    flist = ",".join(f"#{f}" for f in face_ids)
    org = add("CARTESIAN_POINT('',(0.,0.,0.))")
    oz = add("DIRECTION('',(0.,0.,1.))")
    ox = add("DIRECTION('',(1.,0.,0.))")
    place = add(f"AXIS2_PLACEMENT_3D('',#{org},#{oz},#{ox})")

    if m.closed:
        shell = add(f"CLOSED_SHELL('',({flist}))")
        body = add(f"MANIFOLD_SOLID_BREP('{name}',#{shell})")
        rep = add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('{name}',"
                  f"(#{place},#{body}),#{ctx})")
    else:
        shell = add(f"OPEN_SHELL('',({flist}))")
        body = add(f"SHELL_BASED_SURFACE_MODEL('{name}',(#{shell}))")
        rep = add(f"MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{name}',"
                  f"(#{place},#{body}),#{ctx})")

    # --- product structure
    appctx = add("APPLICATION_CONTEXT('core data for automotive mechanical "
                 "design processes')")
    add(f"APPLICATION_PROTOCOL_DEFINITION('international standard',"
        f"'automotive_design',2000,#{appctx})")
    pctx = add(f"PRODUCT_CONTEXT('',#{appctx},'mechanical')")
    prod = add(f"PRODUCT('{name}','{name}','',(#{pctx}))")
    add(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{prod}))")
    pdf = add(f"PRODUCT_DEFINITION_FORMATION('','',#{prod})")
    pdctx = add(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{appctx},"
                "'design')")
    pdef = add(f"PRODUCT_DEFINITION('design','',#{pdf},#{pdctx})")
    pds = add(f"PRODUCT_DEFINITION_SHAPE('','',#{pdef})")
    add(f"SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})")

    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("ISO-10303-21;\nHEADER;\n")
        f.write("FILE_DESCRIPTION(('Converted from STL by stl2step'),"
                "'2;1');\n")
        f.write(f"FILE_NAME('{_q(os.path.basename(path))}','{now}',('')"
                ",(''),'stl2step','stl2step','');\n")
        f.write("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'"
                "));\nENDSEC;\nDATA;\n")
        for i, body_s in enumerate(S.lines, 1):
            f.write(f"#{i}={body_s};\n")
        f.write("ENDSEC;\nEND-ISO-10303-21;\n")
    return len(face_ids), len(edges)


# --------------------------------------------------------------------------
# top-level conversion
# --------------------------------------------------------------------------

def convert(src, dst, merge=True, scale=1.0, weld_rel=1e-5,
            angle_tol=0.1, log=print):
    log(f"Reading {src} ...")
    tris = read_stl(src)
    log(f"  {len(tris):,} triangles")
    m = build_mesh(tris, weld_rel=weld_rel, scale=scale)
    log(f"  {len(m.verts):,} unique vertices, {len(m.faces):,} valid "
        f"triangles, {'watertight' if m.closed else 'NOT watertight'}")
    for w in m.warnings:
        log("  ! " + w)

    if merge and not m.nonmanifold:
        log("Merging coplanar triangles ...")
        geo, merged = merge_planar(m, ang_tol_deg=angle_tol)
        log(f"  {len(m.faces):,} triangles -> {len(geo):,} faces "
            f"({merged:,} triangles merged into planar faces)")
    else:
        geo = [(_tri_normal(m.verts, f), list(f), []) for f in m.faces]
        if merge:
            log("  ! Skipping merge because of non-manifold edges.")

    log("Writing STEP ...")
    nf, ne = write_step(m, geo, dst,
                        name=os.path.splitext(os.path.basename(src))[0])
    size = os.path.getsize(dst) / 1e6
    log(f"Done: {nf:,} faces, {ne:,} edges, {size:.1f} MB -> {dst}")
    return m.closed


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------

def run_gui():
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    root = tk.Tk()
    root.title("STL -> STEP converter")
    root.geometry("640x520")

    src_v = tk.StringVar()
    dst_v = tk.StringVar()
    merge_v = tk.BooleanVar(value=True)
    scale_v = tk.StringVar(value="1.0")
    ang_v = tk.StringVar(value="0.1")

    frm = ttk.Frame(root, padding=12)
    frm.pack(fill="both", expand=True)
    frm.columnconfigure(1, weight=1)

    def pick_src():
        p = filedialog.askopenfilename(
            filetypes=[("STL files", "*.stl *.STL"), ("All files", "*.*")])
        if p:
            src_v.set(p)
            dst_v.set(os.path.splitext(p)[0] + ".step")

    def pick_dst():
        p = filedialog.asksaveasfilename(
            defaultextension=".step",
            filetypes=[("STEP files", "*.step *.stp")])
        if p:
            dst_v.set(p)

    ttk.Label(frm, text="STL file").grid(row=0, column=0, sticky="w", pady=4)
    ttk.Entry(frm, textvariable=src_v).grid(row=0, column=1, sticky="ew",
                                            padx=6)
    ttk.Button(frm, text="Browse...", command=pick_src).grid(row=0, column=2)
    ttk.Label(frm, text="STEP output").grid(row=1, column=0, sticky="w",
                                            pady=4)
    ttk.Entry(frm, textvariable=dst_v).grid(row=1, column=1, sticky="ew",
                                            padx=6)
    ttk.Button(frm, text="Save as...", command=pick_dst).grid(row=1, column=2)

    ttk.Checkbutton(frm, text="Merge coplanar triangles into planar faces "
                    "(recommended)", variable=merge_v).grid(
        row=2, column=0, columnspan=3, sticky="w", pady=4)
    ttk.Label(frm, text="Scale factor").grid(row=3, column=0, sticky="w")
    ttk.Entry(frm, textvariable=scale_v, width=10).grid(row=3, column=1,
                                                        sticky="w", padx=6)
    ttk.Label(frm, text="(STL in inches? use 25.4. Output is millimetres.)"
              ).grid(row=3, column=1, sticky="e")
    ttk.Label(frm, text="Planar tolerance (deg)").grid(row=4, column=0,
                                                       sticky="w", pady=4)
    ttk.Entry(frm, textvariable=ang_v, width=10).grid(row=4, column=1,
                                                      sticky="w", padx=6)

    btn = ttk.Button(frm, text="Convert")
    btn.grid(row=5, column=0, columnspan=3, pady=8, sticky="ew")
    logbox = scrolledtext.ScrolledText(frm, height=16, state="disabled")
    logbox.grid(row=6, column=0, columnspan=3, sticky="nsew")
    frm.rowconfigure(6, weight=1)

    def log(msg):
        def _do():
            logbox.configure(state="normal")
            logbox.insert("end", msg + "\n")
            logbox.see("end")
            logbox.configure(state="disabled")
        root.after(0, _do)

    def work():
        try:
            convert(src_v.get(), dst_v.get(), merge=merge_v.get(),
                    scale=float(scale_v.get()),
                    angle_tol=float(ang_v.get()), log=log)
        except Exception as e:  # noqa
            log(f"ERROR: {e}")
            root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            root.after(0, lambda: btn.configure(state="normal"))

    def go():
        if not src_v.get() or not dst_v.get():
            messagebox.showinfo("Missing file", "Choose an STL file first.")
            return
        btn.configure(state="disabled")
        threading.Thread(target=work, daemon=True).start()

    btn.configure(command=go)
    log("Note: STL has only triangles. Flat areas become real planar faces;\n"
        "curved surfaces stay faceted. This gives you a solid body that\n"
        "SolidWorks can import, but not parametric features.")
    root.mainloop()


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Convert STL to STEP (AP214).")
    ap.add_argument("input", nargs="?", help="input .stl (omit for GUI)")
    ap.add_argument("-o", "--output", help="output .step path")
    ap.add_argument("--no-merge", action="store_true",
                    help="keep one face per triangle")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply coordinates (25.4 for inch STL)")
    ap.add_argument("--angle-tol", type=float, default=0.1,
                    help="coplanarity angle tolerance in degrees")
    args = ap.parse_args()
    if not args.input:
        run_gui()
        return
    out = args.output or os.path.splitext(args.input)[0] + ".step"
    convert(args.input, out, merge=not args.no_merge, scale=args.scale,
            angle_tol=args.angle_tol)


if __name__ == "__main__":
    main()
