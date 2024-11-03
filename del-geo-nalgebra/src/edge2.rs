//! functions for 2D edge (line segment)

use num_traits::AsPrimitive;

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
    let area1 = del_geo_core::tri2::area(po_s0.as_ref(), po_e0.as_ref(), po_s1.as_ref());
    let area2 = del_geo_core::tri2::area(po_s0.as_ref(), po_e0.as_ref(), po_e1.as_ref());
    let area3 = del_geo_core::tri2::area(po_s1.as_ref(), po_e1.as_ref(), po_s0.as_ref());
    let area4 = del_geo_core::tri2::area(po_s1.as_ref(), po_e1.as_ref(), po_e0.as_ref());
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
    if del_geo_core::edge2::intersection_edge2(
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

pub fn barycentric<T>(
    ps: &nalgebra::Vector2<T>,
    pe: &nalgebra::Vector2<T>,
    q: &nalgebra::Vector2<T>,
) -> T
where
    T: nalgebra::RealField + Copy,
{
    let v = pe - ps;
    let t = v.dot(q);
    let s = v.dot(&v);
    t / s.sqrt()
}
