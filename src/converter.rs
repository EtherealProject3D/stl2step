//! STL -> STEP (AP214) conversion core. Standard library only.
//!
//! Pipeline: read STL -> weld vertices / clean -> (optional) merge coplanar
//! triangles into planar faces with holes -> write STEP B-rep.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

pub type V3 = [f64; 3];

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
fn unit(a: V3) -> V3 {
    let n = norm(a);
    if n > 0.0 {
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        [0.0, 0.0, 1.0]
    }
}

pub struct Options {
    /// merge coplanar neighbouring triangles into planar faces
    pub merge: bool,
    /// multiply all coordinates by this factor (25.4 for inch STLs)
    pub scale: f64,
    /// vertex weld tolerance relative to the bounding-box diagonal
    pub weld_rel: f64,
    /// coplanarity angle tolerance in degrees
    pub angle_tol: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            merge: true,
            scale: 1.0,
            weld_rel: 1e-5,
            angle_tol: 0.1,
        }
    }
}

// ---------------------------------------------------------------------------
// STL reading
// ---------------------------------------------------------------------------

pub fn read_stl(path: &Path) -> Result<Vec<[V3; 3]>, String> {
    let data = std::fs::read(path).map_err(|e| format!("Cannot read file: {e}"))?;
    if data.len() < 15 {
        return Err("File is too small to be an STL.".into());
    }

    let mut binary = false;
    if data.len() >= 84 {
        let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
        if 84usize.saturating_add(n.saturating_mul(50)) == data.len() {
            binary = true;
        } else if !data[..5].eq_ignore_ascii_case(b"solid") {
            binary = true; // binary with a wrong/truncated count
        }
    }

    let mut tris: Vec<[V3; 3]> = Vec::new();
    if binary {
        let declared = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
        let n = declared.min((data.len() - 84) / 50);
        tris.reserve(n);
        let f = |off: usize| -> f64 {
            f32::from_le_bytes([data[off], data[off + 1], data[off + 2], data[off + 3]]) as f64
        };
        for i in 0..n {
            let o = 84 + i * 50 + 12; // skip the 12-byte normal
            let mut t = [[0.0f64; 3]; 3];
            for v in 0..3 {
                for k in 0..3 {
                    t[v][k] = f(o + (v * 3 + k) * 4);
                }
            }
            tris.push(t);
        }
    } else {
        let text = String::from_utf8_lossy(&data);
        let mut toks = text.split_whitespace();
        let mut vs: Vec<V3> = Vec::new();
        while let Some(t) = toks.next() {
            if t.eq_ignore_ascii_case("vertex") {
                let mut v = [0.0f64; 3];
                for k in 0..3 {
                    let s = toks.next().ok_or("Truncated vertex in ASCII STL.")?;
                    v[k] = s
                        .parse::<f64>()
                        .map_err(|_| format!("Bad number '{s}' in ASCII STL."))?;
                }
                vs.push(v);
            }
        }
        if vs.len() % 3 != 0 {
            return Err("ASCII STL vertex count is not divisible by 3.".into());
        }
        for c in vs.chunks_exact(3) {
            tris.push([c[0], c[1], c[2]]);
        }
    }
    if tris.is_empty() {
        return Err("No triangles found in file.".into());
    }
    Ok(tris)
}

// ---------------------------------------------------------------------------
// Mesh cleanup
// ---------------------------------------------------------------------------

pub struct Mesh {
    pub verts: Vec<V3>,
    pub faces: Vec<[usize; 3]>,
    pub closed: bool,
    /// true if any edge is shared by 3+ triangles
    pub nonmanifold: bool,
    pub open_edges: usize,
    pub nm_edges: usize,
    pub diag: f64,
    pub warnings: Vec<String>,
}

fn canon(a: usize, b: usize, c: usize) -> [usize; 3] {
    if a <= b && a <= c {
        [a, b, c]
    } else if b <= a && b <= c {
        [b, c, a]
    } else {
        [c, a, b]
    }
}

#[inline]
fn ukey(u: usize, v: usize) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

/// Flip triangles so that neighbouring triangles (across clean 2-triangle
/// edges) agree on winding. Returns the number of triangles flipped.
fn orient_consistent(faces: &mut [[usize; 3]]) -> usize {
    let n = faces.len();
    let mut emap: HashMap<(usize, usize), Vec<usize>> = HashMap::with_capacity(n * 2);
    for (t, f) in faces.iter().enumerate() {
        for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            emap.entry(ukey(u, v)).or_default().push(t);
        }
    }
    let mut seen = vec![false; n];
    let mut flips = 0usize;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        seen[s] = true;
        let mut stack = vec![s];
        while let Some(t) = stack.pop() {
            let f = faces[t];
            for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                let lst = &emap[&ukey(u, v)];
                if lst.len() != 2 {
                    continue;
                }
                let nb = if lst[0] == t { lst[1] } else { lst[0] };
                if seen[nb] {
                    continue;
                }
                seen[nb] = true;
                let g = faces[nb];
                let has_uv = (g[0] == u && g[1] == v)
                    || (g[1] == u && g[2] == v)
                    || (g[2] == u && g[0] == v);
                if has_uv {
                    faces[nb] = [g[0], g[2], g[1]];
                    flips += 1;
                }
                stack.push(nb);
            }
        }
    }
    flips
}

pub fn build_mesh(tris: &[[V3; 3]], weld_rel: f64, scale: f64) -> Result<Mesh, String> {
    let mut mn = [f64::INFINITY; 3];
    let mut mx = [f64::NEG_INFINITY; 3];
    for t in tris {
        for p in t.iter() {
            for k in 0..3 {
                let c = p[k] * scale;
                if !c.is_finite() {
                    return Err("STL contains NaN or infinite coordinates.".into());
                }
                mn[k] = mn[k].min(c);
                mx[k] = mx[k].max(c);
            }
        }
    }
    let d = norm(sub(mx, mn));
    let diag = if d > 0.0 { d } else { 1.0 };
    let tol = diag * weld_rel;
    let inv = 1.0 / tol;

    let mut verts: Vec<V3> = Vec::new();
    let mut vid: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut raw: Vec<[usize; 3]> = Vec::with_capacity(tris.len());
    let mut degenerate = 0usize;
    let min_cross2 = 4.0 * (1e-12 * diag * diag) * (1e-12 * diag * diag);

    for t in tris {
        let mut idx = [0usize; 3];
        for (j, p) in t.iter().enumerate() {
            let q = [p[0] * scale, p[1] * scale, p[2] * scale];
            let key = (
                (q[0] * inv).round() as i64,
                (q[1] * inv).round() as i64,
                (q[2] * inv).round() as i64,
            );
            let i = *vid.entry(key).or_insert_with(|| {
                verts.push(q);
                verts.len() - 1
            });
            idx[j] = i;
        }
        let (a, b, c) = (idx[0], idx[1], idx[2]);
        if a == b || b == c || a == c {
            degenerate += 1;
            continue;
        }
        let cr = cross(sub(verts[b], verts[a]), sub(verts[c], verts[a]));
        if dot(cr, cr) < min_cross2 {
            degenerate += 1;
            continue;
        }
        raw.push([a, b, c]);
    }

    // Duplicates: same-orientation copies collapse to one; back-to-back pairs
    // (same vertices, opposite winding = internal double wall) cancel out.
    let mut gidx: HashMap<[usize; 3], usize> = HashMap::new();
    let mut groups: Vec<(usize, usize, [usize; 3])> = Vec::new();
    for f in &raw {
        let mut key = *f;
        key.sort_unstable();
        let rot = canon(f[0], f[1], f[2]);
        match gidx.get(&key) {
            None => {
                gidx.insert(key, groups.len());
                groups.push((1, 0, rot));
            }
            Some(&gi) => {
                if groups[gi].2 == rot {
                    groups[gi].0 += 1;
                } else {
                    groups[gi].1 += 1;
                }
            }
        }
    }
    drop(raw);
    let mut faces: Vec<[usize; 3]> = Vec::with_capacity(groups.len());
    let mut removed = 0usize;
    for &(npos, nneg, rf) in &groups {
        let cancel = npos.min(nneg);
        let p = npos - cancel;
        let q = nneg - cancel;
        let keep = if p > 0 || q > 0 { 1 } else { 0 };
        removed += npos + nneg - keep;
        if p > 0 {
            faces.push(rf);
        } else if q > 0 {
            faces.push([rf[0], rf[2], rf[1]]);
        }
    }

    let mut warnings = Vec::new();
    if degenerate > 0 {
        warnings.push(format!("Removed {degenerate} degenerate triangle(s)."));
    }
    if removed > 0 {
        warnings.push(format!(
            "Removed {removed} duplicate / back-to-back triangle(s)."
        ));
    }
    if faces.is_empty() {
        return Err("No valid triangles left after cleanup.".into());
    }

    let flips = orient_consistent(&mut faces);
    if flips > 0 {
        warnings.push(format!("Fixed inconsistent winding on {flips} triangle(s)."));
    }

    // topology
    let mut ucount: HashMap<(usize, usize), u32> = HashMap::new();
    let mut dcount: HashMap<(usize, usize), u32> = HashMap::new();
    for f in &faces {
        for (a, b) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            *ucount.entry(ukey(a, b)).or_insert(0) += 1;
            *dcount.entry((a, b)).or_insert(0) += 1;
        }
    }
    let open_edges = ucount.values().filter(|&&c| c == 1).count();
    let nm_edges = ucount.values().filter(|&&c| c > 2).count();
    let nonmanifold = nm_edges > 0;
    let closed = open_edges == 0
        && nm_edges == 0
        && dcount
            .iter()
            .all(|(&(a, b), &c)| c == 1 && dcount.contains_key(&(b, a)));

    if closed {
        let mut vol = 0.0;
        for f in &faces {
            vol += dot(verts[f[0]], cross(verts[f[1]], verts[f[2]]));
        }
        if vol < 0.0 {
            for f in faces.iter_mut() {
                f.swap(1, 2);
            }
            warnings.push("Mesh was inside-out; flipped it.".into());
        }
    } else {
        warnings.push(format!(
            "Mesh is not a clean closed solid ({open_edges} open edge(s), {nm_edges} \
             non-manifold edge(s)). Output will be a surface model; in SolidWorks use \
             Knit Surface with 'Try to form solid'. Fixing the STL first gives a true solid."
        ));
    }

    Ok(Mesh {
        verts,
        faces,
        closed,
        nonmanifold,
        open_edges,
        nm_edges,
        diag,
        warnings,
    })
}

// ---------------------------------------------------------------------------
// Planar face merging
// ---------------------------------------------------------------------------

pub struct FaceGeo {
    pub normal: V3,
    pub outer: Vec<usize>,
    pub holes: Vec<Vec<usize>>,
}

fn tri_normal(verts: &[V3], f: &[usize; 3]) -> V3 {
    let (a, b, c) = (verts[f[0]], verts[f[1]], verts[f[2]]);
    unit(cross(sub(b, a), sub(c, a)))
}

fn single_faces(m: &Mesh) -> Vec<FaceGeo> {
    m.faces
        .iter()
        .map(|f| FaceGeo {
            normal: tri_normal(&m.verts, f),
            outer: vec![f[0], f[1], f[2]],
            holes: Vec::new(),
        })
        .collect()
}

/// Returns the face list and how many triangles were merged into larger faces.
pub fn merge_planar(m: &Mesh, ang_tol_deg: f64, dist_rel: f64) -> (Vec<FaceGeo>, usize) {
    let verts = &m.verts;
    let faces = &m.faces;
    let n = faces.len();
    let normals: Vec<V3> = faces.iter().map(|f| tri_normal(verts, f)).collect();
    let cos_tol = ang_tol_deg.to_radians().cos();
    let dist_tol = dist_rel * m.diag;

    let mut ucount: HashMap<(usize, usize), u32> = HashMap::with_capacity(n * 2);
    for f in faces.iter() {
        for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            *ucount.entry(ukey(u, v)).or_insert(0) += 1;
        }
    }
    // Only cross clean edges shared by exactly two triangles.
    let mut emap: HashMap<(usize, usize), usize> = HashMap::with_capacity(n * 3);
    for (t, f) in faces.iter().enumerate() {
        for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
            if ucount[&ukey(u, v)] == 2 {
                emap.insert((u, v), t);
            }
        }
    }

    let mut visited = vec![false; n];
    let mut out: Vec<FaceGeo> = Vec::new();
    let mut merged = 0usize;

    for seed in 0..n {
        if visited[seed] {
            continue;
        }
        let n0 = normals[seed];
        let d0 = dot(n0, verts[faces[seed][0]]);
        visited[seed] = true;
        let mut region = vec![seed];
        let mut stack = vec![seed];
        while let Some(t) = stack.pop() {
            let f = faces[t];
            for (u, v) in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                if let Some(&nb) = emap.get(&(v, u)) {
                    if visited[nb] {
                        continue;
                    }
                    if dot(normals[nb], n0) < cos_tol {
                        continue;
                    }
                    if faces[nb]
                        .iter()
                        .any(|&i| (dot(n0, verts[i]) - d0).abs() > dist_tol)
                    {
                        continue;
                    }
                    visited[nb] = true;
                    region.push(nb);
                    stack.push(nb);
                }
            }
        }

        let loops = if region.len() > 1 {
            region_loops(&region, faces, verts, n0)
        } else {
            None
        };
        match loops {
            Some(mut ls) => {
                merged += region.len();
                let outer = ls.remove(0);
                out.push(FaceGeo {
                    normal: n0,
                    outer,
                    holes: ls,
                });
            }
            None => {
                for &t in &region {
                    let f = faces[t];
                    out.push(FaceGeo {
                        normal: normals[t],
                        outer: vec![f[0], f[1], f[2]],
                        holes: Vec::new(),
                    });
                }
            }
        }
    }
    (out, merged)
}

/// Boundary loops of a triangle region. First loop is the outer one (CCW about
/// n0), the rest are holes (CW). Returns None if the region is not clean
/// (e.g. pinch points); the caller then falls back to individual triangles.
fn region_loops(
    region: &[usize],
    faces: &[[usize; 3]],
    verts: &[V3],
    n0: V3,
) -> Option<Vec<Vec<usize>>> {
    let mut dset: HashSet<(usize, usize)> = HashSet::new();
    for &t in region {
        let f = faces[t];
        dset.insert((f[0], f[1]));
        dset.insert((f[1], f[2]));
        dset.insert((f[2], f[0]));
    }
    if dset.len() != region.len() * 3 {
        return None; // overlapping triangles; fall back
    }
    let mut nxt: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(a, b) in &dset {
        if !dset.contains(&(b, a)) {
            nxt.entry(a).or_default().push(b);
        }
    }
    let mut next1: HashMap<usize, usize> = HashMap::new();
    for (a, v) in &nxt {
        if v.len() != 1 {
            return None;
        }
        next1.insert(*a, v[0]);
    }
    let mut keys: Vec<usize> = next1.keys().copied().collect();
    keys.sort_unstable();

    let mut seen: HashSet<usize> = HashSet::new();
    let mut loops: Vec<Vec<usize>> = Vec::new();
    for &s in &keys {
        if seen.contains(&s) {
            continue;
        }
        let mut lp = vec![s];
        seen.insert(s);
        let mut cur = next1[&s];
        while cur != s {
            if seen.contains(&cur) || !next1.contains_key(&cur) {
                return None;
            }
            lp.push(cur);
            seen.insert(cur);
            cur = next1[&cur];
        }
        loops.push(lp);
    }

    let areas: Vec<f64> = loops
        .iter()
        .map(|lp| {
            let mut tot = 0.0;
            for i in 0..lp.len() {
                let p = verts[lp[i]];
                let q = verts[lp[(i + 1) % lp.len()]];
                tot += dot(cross(p, q), n0);
            }
            0.5 * tot
        })
        .collect();
    let pos: Vec<usize> = (0..loops.len()).filter(|&i| areas[i] > 0.0).collect();
    if pos.len() != 1 {
        return None;
    }
    let o = pos[0];
    let mut res = vec![loops[o].clone()];
    for (i, l) in loops.iter().enumerate() {
        if i != o {
            res.push(l.clone());
        }
    }
    Some(res)
}

// ---------------------------------------------------------------------------
// STEP writer
// ---------------------------------------------------------------------------

/// Format a float as a STEP REAL (must contain a decimal point).
fn rr(x: f64) -> String {
    let s = format!("{:?}", x);
    if let Some(pos) = s.find(|c: char| c == 'e' || c == 'E') {
        let mant = &s[..pos];
        let exp = &s[pos + 1..];
        let mant = if mant.contains('.') {
            mant.to_string()
        } else {
            format!("{mant}.")
        };
        format!("{mant}E{exp}")
    } else if s.contains('.') {
        s
    } else {
        format!("{s}.")
    }
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn iso_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    // civil-from-days (H. Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let mut y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    if mo <= 2 {
        y += 1;
    }
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}",
        y,
        mo,
        d,
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
}

struct Builder<'a> {
    lines: Vec<String>,
    verts: &'a [V3],
    pt: Vec<usize>,
    vx: Vec<usize>,
    edges: HashMap<(usize, usize), usize>,
}

impl<'a> Builder<'a> {
    fn add(&mut self, s: String) -> usize {
        self.lines.push(s);
        self.lines.len() // STEP ids are 1-based
    }

    fn edge_curve(&mut self, a: usize, b: usize) -> usize {
        let key = if a < b { (a, b) } else { (b, a) };
        if let Some(&e) = self.edges.get(&key) {
            return e;
        }
        let (lo, hi) = key;
        let d = sub(self.verts[hi], self.verts[lo]);
        let len = norm(d);
        let u = unit(d);
        let di = self.add(format!(
            "DIRECTION('',({},{},{}))",
            rr(u[0]),
            rr(u[1]),
            rr(u[2])
        ));
        let ve = self.add(format!("VECTOR('',#{},{})", di, rr(len)));
        let ln = self.add(format!("LINE('',#{},#{})", self.pt[lo], ve));
        let ec = self.add(format!(
            "EDGE_CURVE('',#{},#{},#{},.T.)",
            self.vx[lo], self.vx[hi], ln
        ));
        self.edges.insert(key, ec);
        ec
    }

    fn loop_id(&mut self, lp: &[usize]) -> usize {
        let k = lp.len();
        let mut oes: Vec<String> = Vec::with_capacity(k);
        for i in 0..k {
            let (a, b) = (lp[i], lp[(i + 1) % k]);
            let ec = self.edge_curve(a, b);
            let oe = self.add(format!(
                "ORIENTED_EDGE('',*,*,#{},{})",
                ec,
                if a < b { ".T." } else { ".F." }
            ));
            oes.push(format!("#{oe}"));
        }
        self.add(format!("EDGE_LOOP('',({}))", oes.join(",")))
    }
}

pub fn write_step(
    m: &Mesh,
    geo: &[FaceGeo],
    path: &Path,
    name: &str,
) -> Result<(usize, usize), String> {
    let name = sanitize(name);
    let mut b = Builder {
        lines: Vec::new(),
        verts: &m.verts,
        pt: Vec::with_capacity(m.verts.len()),
        vx: Vec::with_capacity(m.verts.len()),
        edges: HashMap::new(),
    };

    // units & context
    let lu = b.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))".into());
    let au = b.add("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))".into());
    let su = b.add("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())".into());
    let unc = b.add(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#{lu},\
         'distance_accuracy_value','confusion accuracy')"
    ));
    let ctx = b.add(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3) \
         GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{unc})) \
         GLOBAL_UNIT_ASSIGNED_CONTEXT((#{lu},#{au},#{su})) \
         REPRESENTATION_CONTEXT('Context3D','3D Context with UNIT and UNCERTAINTY'))"
    ));

    // vertices
    for i in 0..m.verts.len() {
        let p = m.verts[i];
        let c = b.add(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            rr(p[0]),
            rr(p[1]),
            rr(p[2])
        ));
        b.pt.push(c);
        let v = b.add(format!("VERTEX_POINT('',#{c})"));
        b.vx.push(v);
    }

    // faces
    let mut face_ids: Vec<usize> = Vec::with_capacity(geo.len());
    for g in geo {
        let p = m.verts[g.outer[0]];
        let o = b.add(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            rr(p[0]),
            rr(p[1]),
            rr(p[2])
        ));
        let nrm = g.normal;
        let mut axis = 0usize;
        for k in 1..3 {
            if nrm[k].abs() < nrm[axis].abs() {
                axis = k;
            }
        }
        let mut helper = [0.0f64; 3];
        helper[axis] = 1.0;
        let rf = unit(cross(nrm, helper));
        let zd = b.add(format!(
            "DIRECTION('',({},{},{}))",
            rr(nrm[0]),
            rr(nrm[1]),
            rr(nrm[2])
        ));
        let xd = b.add(format!(
            "DIRECTION('',({},{},{}))",
            rr(rf[0]),
            rr(rf[1]),
            rr(rf[2])
        ));
        let ax = b.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{zd},#{xd})"));
        let pl = b.add(format!("PLANE('',#{ax})"));

        let outer_loop = b.loop_id(&g.outer);
        let mut bounds = vec![b.add(format!("FACE_OUTER_BOUND('',#{outer_loop},.T.)"))];
        for h in &g.holes {
            let hl = b.loop_id(h);
            bounds.push(b.add(format!("FACE_BOUND('',#{hl},.T.)")));
        }
        let bl: Vec<String> = bounds.iter().map(|x| format!("#{x}")).collect();
        let f = b.add(format!("ADVANCED_FACE('',({}),#{pl},.T.)", bl.join(",")));
        face_ids.push(f);
    }
    let flist: Vec<String> = face_ids.iter().map(|x| format!("#{x}")).collect();
    let flist = flist.join(",");

    let org = b.add("CARTESIAN_POINT('',(0.,0.,0.))".into());
    let oz = b.add("DIRECTION('',(0.,0.,1.))".into());
    let ox = b.add("DIRECTION('',(1.,0.,0.))".into());
    let place = b.add(format!("AXIS2_PLACEMENT_3D('',#{org},#{oz},#{ox})"));

    let rep;
    if m.closed {
        let shell = b.add(format!("CLOSED_SHELL('',({flist}))"));
        let body = b.add(format!("MANIFOLD_SOLID_BREP('{name}',#{shell})"));
        rep = b.add(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{name}',(#{place},#{body}),#{ctx})"
        ));
    } else {
        let shell = b.add(format!("OPEN_SHELL('',({flist}))"));
        let body = b.add(format!("SHELL_BASED_SURFACE_MODEL('{name}',(#{shell}))"));
        rep = b.add(format!(
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{name}',(#{place},#{body}),#{ctx})"
        ));
    }

    // product structure
    let appctx = b.add(
        "APPLICATION_CONTEXT('core data for automotive mechanical design processes')".into(),
    );
    b.add(format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{appctx})"
    ));
    let pctx = b.add(format!("PRODUCT_CONTEXT('',#{appctx},'mechanical')"));
    let prod = b.add(format!("PRODUCT('{name}','{name}','',(#{pctx}))"));
    b.add(format!("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{prod}))"));
    let pdf = b.add(format!("PRODUCT_DEFINITION_FORMATION('','',#{prod})"));
    let pdctx = b.add(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{appctx},'design')"
    ));
    let pdef = b.add(format!("PRODUCT_DEFINITION('design','',#{pdf},#{pdctx})"));
    let pds = b.add(format!("PRODUCT_DEFINITION_SHAPE('','',#{pdef})"));
    b.add(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));

    // file
    let file_name = path
        .file_name()
        .map(|s| sanitize(&s.to_string_lossy()))
        .unwrap_or_else(|| "out.step".into());
    let mut out = String::with_capacity(b.lines.len() * 60 + 1024);
    out.push_str("ISO-10303-21;\nHEADER;\n");
    out.push_str("FILE_DESCRIPTION(('Converted from STL by stl2step'),'2;1');\n");
    out.push_str(&format!(
        "FILE_NAME('{}','{}',(''),(''),'stl2step','stl2step','');\n",
        file_name,
        iso_now()
    ));
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\nDATA;\n");
    for (i, l) in b.lines.iter().enumerate() {
        out.push_str(&format!("#{}={};\n", i + 1, l));
    }
    out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");

    std::fs::write(path, out.as_bytes()).map_err(|e| format!("Cannot write output: {e}"))?;
    Ok((face_ids.len(), b.edges.len()))
}

// ---------------------------------------------------------------------------
// Top-level
// ---------------------------------------------------------------------------

pub fn convert(
    src: &Path,
    dst: &Path,
    opts: &Options,
    log: &mut dyn FnMut(String),
) -> Result<bool, String> {
    if !(opts.scale.is_finite() && opts.scale > 0.0) {
        return Err("Scale must be a positive number.".into());
    }
    log(format!("Reading {} ...", src.display()));
    let tris = read_stl(src)?;
    log(format!("  {} triangles", tris.len()));

    let m = build_mesh(&tris, opts.weld_rel, opts.scale)?;
    drop(tris);
    log(format!(
        "  {} unique vertices, {} valid triangles, {}",
        m.verts.len(),
        m.faces.len(),
        if m.closed { "watertight" } else { "NOT watertight" }
    ));
    for w in &m.warnings {
        log(format!("  ! {w}"));
    }

    let geo = if opts.merge {
        log("Merging coplanar triangles ...".into());
        let (g, merged) = merge_planar(&m, opts.angle_tol, 1e-4);
        log(format!(
            "  {} triangles -> {} faces ({} triangles merged into planar faces)",
            m.faces.len(),
            g.len(),
            merged
        ));
        g
    } else {
        single_faces(&m)
    };

    if geo.len() > 20_000 {
        log(format!(
            "  ! {} faces is a LOT for SolidWorks and may take very long to open. \
             Decimate the STL (e.g. in MeshLab) first.",
            geo.len()
        ));
    }

    log("Writing STEP ...".into());
    let stem = src
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "part".into());
    let (nf, ne) = write_step(&m, &geo, dst, &stem)?;
    let size = std::fs::metadata(dst).map(|md| md.len()).unwrap_or(0) as f64 / 1e6;
    log(format!(
        "Done: {} faces, {} edges, {:.1} MB -> {}",
        nf,
        ne,
        size,
        dst.display()
    ));
    Ok(m.closed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit cube, 12 triangles, outward-facing CCW winding.
    fn cube() -> Vec<[V3; 3]> {
        vec![
            // +z
            [[0., 0., 1.], [1., 0., 1.], [1., 1., 1.]],
            [[0., 0., 1.], [1., 1., 1.], [0., 1., 1.]],
            // -z
            [[0., 0., 0.], [1., 1., 0.], [1., 0., 0.]],
            [[0., 0., 0.], [0., 1., 0.], [1., 1., 0.]],
            // +x
            [[1., 0., 0.], [1., 1., 0.], [1., 1., 1.]],
            [[1., 0., 0.], [1., 1., 1.], [1., 0., 1.]],
            // -x
            [[0., 0., 0.], [0., 1., 1.], [0., 1., 0.]],
            [[0., 0., 0.], [0., 0., 1.], [0., 1., 1.]],
            // +y
            [[0., 1., 0.], [0., 1., 1.], [1., 1., 1.]],
            [[0., 1., 0.], [1., 1., 1.], [1., 1., 0.]],
            // -y
            [[0., 0., 0.], [1., 0., 0.], [1., 0., 1.]],
            [[0., 0., 0.], [1., 0., 1.], [0., 0., 1.]],
        ]
    }

    #[test]
    fn cube_is_closed_and_welded() {
        let m = build_mesh(&cube(), 1e-5, 1.0).unwrap();
        assert_eq!(m.verts.len(), 8);
        assert_eq!(m.faces.len(), 12);
        assert!(m.closed);
        assert!(!m.nonmanifold);
    }

    #[test]
    fn inverted_cube_gets_flipped() {
        let inv: Vec<[V3; 3]> = cube().iter().map(|t| [t[0], t[2], t[1]]).collect();
        let m = build_mesh(&inv, 1e-5, 1.0).unwrap();
        assert!(m.closed);
        let mut vol = 0.0;
        for f in &m.faces {
            vol += dot(m.verts[f[0]], cross(m.verts[f[1]], m.verts[f[2]]));
        }
        assert!(vol > 0.0);
    }

    #[test]
    fn cube_merges_to_six_faces() {
        let m = build_mesh(&cube(), 1e-5, 1.0).unwrap();
        let (geo, merged) = merge_planar(&m, 0.1, 1e-4);
        assert_eq!(geo.len(), 6);
        assert_eq!(merged, 12);
        assert!(geo.iter().all(|g| g.outer.len() == 4 && g.holes.is_empty()));
    }

    #[test]
    fn mixed_winding_is_repaired_and_merges() {
        let mut t = cube();
        for i in [0usize, 3, 4, 9] {
            let x = t[i];
            t[i] = [x[0], x[2], x[1]];
        }
        let m = build_mesh(&t, 1e-5, 1.0).unwrap();
        assert!(m.closed);
        let (geo, _) = merge_planar(&m, 0.1, 1e-4);
        assert_eq!(geo.len(), 6);
    }

    #[test]
    fn back_to_back_pair_and_duplicate_are_removed() {
        let mut t = cube();
        let dup = t[1];
        t.push(dup); // same-orientation duplicate
        let a = [[0., 0., 0.], [1., 1., 0.], [1., 1., 1.]];
        t.push(a);
        t.push([a[0], a[2], a[1]]); // internal double wall
        let m = build_mesh(&t, 1e-5, 1.0).unwrap();
        assert_eq!(m.faces.len(), 12);
        assert!(m.closed);
    }

    #[test]
    fn fin_creates_nonmanifold_edge_but_cube_still_merges() {
        let mut t = cube();
        t.push([[0., 0., 0.], [1., 0., 0.], [0.5, -1., 0.5]]);
        let m = build_mesh(&t, 1e-5, 1.0).unwrap();
        assert!(!m.closed);
        assert_eq!(m.nm_edges, 1);
        let (geo, _) = merge_planar(&m, 0.1, 1e-4);
        assert_eq!(geo.len(), 7);
    }

    #[test]
    fn open_mesh_is_not_closed() {
        let mut t = cube();
        t.pop();
        let m = build_mesh(&t, 1e-5, 1.0).unwrap();
        assert!(!m.closed);
    }

    #[test]
    fn writes_valid_looking_step() {
        let m = build_mesh(&cube(), 1e-5, 1.0).unwrap();
        let (geo, _) = merge_planar(&m, 0.1, 1e-4);
        let path = std::env::temp_dir().join("stl2step_test_cube.step");
        let (nf, ne) = write_step(&m, &geo, &path, "cube").unwrap();
        assert_eq!(nf, 6);
        assert_eq!(ne, 12);
        let s = std::fs::read_to_string(&path).unwrap();
        assert!(s.starts_with("ISO-10303-21;"));
        assert!(s.contains("MANIFOLD_SOLID_BREP"));
        assert!(s.contains("CLOSED_SHELL"));
        assert!(s.trim_end().ends_with("END-ISO-10303-21;"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn real_formatting_has_decimal_point() {
        assert_eq!(rr(10.0), "10.0");
        assert!(rr(1e-16).contains('.'));
    }
}
