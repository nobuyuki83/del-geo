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

// above: w/o nalgebra
// --------------------------------------------------
// below: use nalgebra

pub fn length_squared<T>(p0: &nalgebra::Vector2<T>, p1: &nalgebra::Vector2<T>) -> T
where
    T: std::ops::Sub<Output = T> + std::ops::Mul<Output = T> + std::ops::Add<Output = T> + Copy,
{
    let x = p0[0] - p1[0];
    let y = p0[1] - p1[1];
    x * x + y * y
}

pub fn intersect_edge2(
    po_s0: &nalgebra::Vector2<f32>,
    po_e0: &nalgebra::Vector2<f32>,
    po_s1: &nalgebra::Vector2<f32>,
    po_e1: &nalgebra::Vector2<f32>,
) -> bool {
    let area1 = crate::tri2::area_(po_s0.as_ref(), po_e0.as_ref(), po_s1.as_ref());
    let area2 = crate::tri2::area_(po_s0.as_ref(), po_e0.as_ref(), po_e1.as_ref());
    let area3 = crate::tri2::area_(po_s1.as_ref(), po_e1.as_ref(), po_s0.as_ref());
    let area4 = crate::tri2::area_(po_s1.as_ref(), po_e1.as_ref(), po_e0.as_ref());
    let a12 = area1 * area2;
    if a12 > 0_f32 {
        return false;
    }
    let a34 = area3 * area4;
    if a34 > 0_f32 {
        return false;
    }
    true
}

pub fn distance_to_edge2<T>(
    po_s0: &nalgebra::Vector2<T>,
    po_e0: &nalgebra::Vector2<T>,
    po_s1: &nalgebra::Vector2<T>,
    po_e1: &nalgebra::Vector2<T>,
) -> T
where
    T: num_traits::Float + nalgebra::RealField + 'static + Copy + std::fmt::Debug,
    f64: AsPrimitive<T>,
{
    if intersection_edge2_(
        po_s0.as_ref(),
        po_e0.as_ref(),
        po_s1.as_ref(),
        po_e1.as_ref(),
    )
    .is_some()
    {
        return (-1_f64).as_();
    }
    let ds1 = crate::edge::distance_to_point(po_s0, po_s1, po_e1);
    let de1 = crate::edge::distance_to_point(po_e0, po_s1, po_e1);
    let ds0 = crate::edge::distance_to_point(po_s1, po_s0, po_e0);
    let de0 = crate::edge::distance_to_point(po_e1, po_s0, po_e0);
    let mut min_dist = ds1;
    min_dist = if de1 < min_dist { de1 } else { min_dist };
    min_dist = if ds0 < min_dist { ds0 } else { min_dist };
    min_dist = if de0 < min_dist { de0 } else { min_dist };
    min_dist
}

pub fn winding_number<T>(
    ps: &nalgebra::Vector2<T>,
    pe: &nalgebra::Vector2<T>,
    po: &nalgebra::Vector2<T>,
) -> T
where
    T: nalgebra::RealField + Copy,
    f64: AsPrimitive<T>,
{
    let p0 = ps - po;
    let p1 = pe - po;
    let y: T = p1[1] * p0[0] - p1[0] * p0[1];
    let x: T = p0[0] * p1[0] + p0[1] * p1[1];
    y.atan2(x) * std::f64::consts::FRAC_1_PI.as_() * 0.5.as_()
}
