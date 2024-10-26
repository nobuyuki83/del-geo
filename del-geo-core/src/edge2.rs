//! methods for 2D edge (line segment)

use num_traits::AsPrimitive;

pub fn length<T>(ps: &[T; 2], pe: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    let dx = ps[0] - pe[0];
    let dy = ps[1] - pe[1];
    (dx * dx + dy * dy).sqrt()
}

pub fn unit_edge_vector<T>(ps: &[T; 2], pe: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    let dx = pe[0] - ps[0];
    let dy = pe[1] - ps[1];
    let linv: T = T::one() / (dx * dx + dy * dy).sqrt();
    [dx * linv, dy * linv]
}

pub fn culling_intersection<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> bool
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let min0x = po_s0[0].min(po_e0[0]);
    let max0x = po_s0[0].max(po_e0[0]);
    let min1x = po_s1[0].min(po_e1[0]);
    let max1x = po_s1[0].max(po_e1[0]);
    let min0y = po_s0[1].min(po_e0[1]);
    let max0y = po_s0[1].max(po_e0[1]);
    let min1y = po_s1[1].min(po_e1[1]);
    let max1y = po_s1[1].max(po_e1[1]);
    let len =
        ((max0x - min0x) + (max0y - min0y) + (max1x - min1x) + (max1y - min1y)) * 0.0001_f64.as_();
    if max1x + len < min0x {
        return false;
    }
    if max0x + len < min1x {
        return false;
    }
    if max1y + len < min0y {
        return false;
    }
    if max0y + len < min1y {
        return false;
    }
    true
}

pub fn intersection_edge2<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> Option<(T, T)>
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let area1 = crate::tri2::area(po_s0, po_e0, po_s1);
    let area2 = crate::tri2::area(po_s0, po_e0, po_e1);
    let area3 = crate::tri2::area(po_s1, po_e1, po_s0);
    let area4 = crate::tri2::area(po_s1, po_e1, po_e0);
    if area1 * area2 > 0_f64.as_() {
        return None;
    }
    if area3 * area4 > 0_f64.as_() {
        return None;
    }
    let r1 = area1 / (area1 - area2);
    let r0 = area3 / (area3 - area4);
    Some((r0, r1))
}

pub fn winding_number<T>(ps: &[T; 2], pe: &[T; 2], po: &[T; 2]) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    let half = T::one() / (T::one() + T::one());
    let p0 = crate::vec2::sub(ps, po);
    let p1 = crate::vec2::sub(pe, po);
    let y: T = p1[1] * p0[0] - p1[0] * p0[1];
    let x: T = p0[0] * p1[0] + p0[1] * p1[1];
    y.atan2(x) * T::FRAC_1_PI() * half
}

/// Find the nearest point on a line segment to the origin(0,0)
/// Returns (k,v), where k is the coeffcient between [0,1], v is the point
pub fn nearest_origin(ps: &[f32; 2], pe: &[f32; 2]) -> (f32, [f32; 2]) {
    let d = crate::vec2::sub(pe, ps);
    let a = crate::vec2::squared_length(&d);
    if a == 0f32 {
        return (0.5, [(ps[0] + pe[0]) * 0.5, (ps[1] + pe[1]) * 0.5]);
    }
    let b = crate::vec2::dot(&d, ps);
    let r0 = (-b / a).clamp(0., 1.);
    (
        r0,
        [
            (1. - r0) * ps[0] + r0 * pe[0],
            (1. - r0) * ps[1] + r0 * pe[1],
        ],
    )
}

#[test]
fn test_nearest_origin() {
    let (_r, pm) = nearest_origin(&[-0.1, 1.0], &[1.0, 1.0]);
    assert!(crate::vec2::length(&crate::vec2::sub(&pm, &[0., 1.])) < 1.0e-5);
}

/// Find the nearest point on a line segment to point p
/// Returns (k,v), where k is the coeffcient, v is the point
pub fn nearest_point2(
    s: &[f32; 2], // source
    e: &[f32; 2], // end
    p: &[f32; 2],
) -> (f32, [f32; 2]) {
    use crate::vec2;
    let (r, p0) = nearest_origin(&vec2::sub(s, p), &vec2::sub(e, p));
    (r, [p0[0] + p[0], p0[1] + p[1]])
}

#[test]
fn test_nearest_point2() {
    let (_r, pm) = nearest_point2(&[-0.1, 1.0], &[1.0, 1.0], &[0.0, 0.3]);
    assert!(crate::vec2::length(&crate::vec2::sub(&pm, &[0., 1.])) < 1.0e-5);
}
