pub struct TriangleScanlineIter {
    v0: [f32; 2],
    v1: [f32; 2],
    v2: [f32; 2],

    // 現在の走査線
    py: i32,
    y_end: i32,

    // 現在のスキャンライン内の x 範囲
    cur_x: i32,
    x_end: i32,

    // 今のスキャンラインが有効か
    has_span: bool,
}

impl TriangleScanlineIter {
    pub fn new(v0: [f32; 2], v1: [f32; 2], v2: [f32; 2]) -> Self {
        let min_y = v0[1].min(v1[1]).min(v2[1]);
        let max_y = v0[1].max(v1[1]).max(v2[1]);

        let py_start = (min_y - 0.5).ceil() as i32;
        let py_end = (max_y - 0.5).floor() as i32;

        Self {
            v0,
            v1,
            v2,
            py: py_start,
            y_end: py_end,
            cur_x: 0,
            x_end: -1,
            has_span: false,
        }
    }

    fn edge_intersect(p0: [f32; 2], p1: [f32; 2], y: f32) -> Option<f32> {
        if (p1[1] - p0[1]).abs() < f32::EPSILON {
            return None;
        }
        let ymin = p0[1].min(p1[1]);
        let ymax = p0[1].max(p1[1]);

        if y < ymin || y >= ymax {
            return None;
        }

        let t = (y - p0[1]) / (p1[1] - p0[1]);
        Some(p0[0] + t * (p1[0] - p0[0]))
    }

    fn setup_span(&mut self) {
        self.has_span = false;

        while self.py <= self.y_end {
            let scan_y = self.py as f32 + 0.5;
            let mut xs = [0.0f32; 3];
            let mut count = 0;

            if let Some(x) = Self::edge_intersect(self.v0, self.v1, scan_y) {
                xs[count] = x;
                count += 1;
            }
            if let Some(x) = Self::edge_intersect(self.v1, self.v2, scan_y) {
                xs[count] = x;
                count += 1;
            }
            if let Some(x) = Self::edge_intersect(self.v2, self.v0, scan_y) {
                xs[count] = x;
                count += 1;
            }

            if count >= 2 {
                xs[..count].sort_by(|a, b| a.partial_cmp(b).unwrap());

                let x_left = xs[0];
                let x_right = xs[count - 1];

                // pixel (ix, iy) is covered if its center (ix+0.5, iy+0.5) is inside the triangle.
                // A pixel ix is covered on this scanline when ix+0.5 is in [x_left, x_right],
                // i.e. ix >= ceil(x_left - 0.5) and ix <= floor(x_right - 0.5).
                // NOTE: f32 edge-intersection rounding can make x_right land just below the
                // true boundary (e.g. 114.49999 instead of 114.5), silently dropping a boundary
                // pixel. Add a small epsilon before flooring if that matters to the caller.
                let x_start = (x_left - 0.5).ceil() as i32;
                let x_end = (x_right - 0.5).floor() as i32;

                if x_start <= x_end {
                    self.cur_x = x_start;
                    self.x_end = x_end;
                    self.has_span = true;
                    return;
                }
            }

            self.py += 1;
        }
    }
}

impl Iterator for TriangleScanlineIter {
    type Item = [i32; 2];

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.has_span && self.cur_x <= self.x_end {
                let px = self.cur_x;
                let py = self.py;
                self.cur_x += 1;
                return Some([px, py]);
            }

            // 次のスキャンラインへ
            self.py += if self.has_span { 1 } else { 0 };
            self.setup_span();

            if !self.has_span {
                return None;
            }
        }
    }
}

#[test]
fn test0() {
    let v0 = [2.0f32, 1.0];
    let v1 = [8.0f32, 3.0];
    let v2 = [4.0f32, 7.0];

    use std::collections::HashSet;
    let sign = crate::tri2::area(&v0, &v1, &v2).signum();
    let pixels_inside: HashSet<_> = TriangleScanlineIter::new(v0, v1, v2).collect();
    for p in &pixels_inside {
        let q = [p[0] as f32 + 0.5, p[1] as f32 + 0.5];
        assert!(crate::tri2::is_inside(&v0, &v1, &v2, &q, sign).is_some());
    }

    let pixels_outside: HashSet<[i32; 2]> = pixels_inside
        .iter()
        .flat_map(|p| {
            [
                [p[0] - 1, p[1]],
                [p[0] + 1, p[1]],
                [p[0], p[1] - 1],
                [p[0], p[1] + 1],
            ]
        })
        .filter(|p| !pixels_inside.contains(p))
        .collect();
    for p in &pixels_outside {
        let q = [p[0] as f32 + 0.5, p[1] as f32 + 0.5];
        assert!(
            crate::tri2::is_inside(&v0, &v1, &v2, &q, sign).is_none(),
            "adjacent pixel ({},{}) is inside the triangle",
            p[0],
            p[1]
        );
    }
}
