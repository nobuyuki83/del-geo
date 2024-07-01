//! methods for 2D edge (line segment)

use num_traits::AsPrimitive;

pub fn culling_intersection_<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> bool
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let min0x = if po_s0[0] < po_e0[0] {
        po_s0[0]
    } else {
        po_e0[0]
    };
    let max0x = if po_s0[0] > po_e0[0] {
        po_s0[0]
    } else {
        po_e0[0]
    };
    let max1x = if po_s1[0] > po_e1[0] {
        po_s1[0]
    } else {
        po_e1[0]
    };
    let min1x = if po_s1[0] < po_e1[0] {
        po_s1[0]
    } else {
        po_e1[0]
    };
    let min0y = if po_s0[1] < po_e0[1] {
        po_s0[1]
    } else {
        po_e0[1]
    };
    let max0y = if po_s0[1] > po_e0[1] {
        po_s0[1]
    } else {
        po_e0[1]
    };
    let max1y = if po_s1[1] > po_e1[1] {
        po_s1[1]
    } else {
        po_e1[1]
    };
    let min1y = if po_s1[1] < po_e1[1] {
        po_s1[1]
    } else {
        po_e1[1]
    };
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

pub fn intersection_edge2_<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> Option<(T, T)>
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let area1 = crate::tri2::area_(po_s0, po_e0, po_s1);
    let area2 = crate::tri2::area_(po_s0, po_e0, po_e1);
    let area3 = crate::tri2::area_(po_s1, po_e1, po_s0);
    let area4 = crate::tri2::area_(po_s1, po_e1, po_e0);
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

pub fn winding_number_<T>(ps: &[T; 2], pe: &[T; 2], po: &[T; 2]) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    let half = T::one() / (T::one() + T::one());
    let p0 = crate::vec2::sub_(ps, po);
    let p1 = crate::vec2::sub_(pe, po);
    let y: T = p1[1] * p0[0] - p1[0] * p0[1];
    let x: T = p0[0] * p1[0] + p0[1] * p1[1];
    y.atan2(x) * T::FRAC_1_PI() * half
}
