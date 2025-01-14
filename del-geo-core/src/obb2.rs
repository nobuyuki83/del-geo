//! 2D Oriented Bounding Box
//! data structure `&[Real;6]`
//! first 2 Reals are for center
//! next 2 Reals are for half of major axis direction
//! next 2 Reals are for half of minar axis direction

use crate::vec2::Vec2;

pub trait OBB2Trait<T>
where
    Self: Sized,
{
    fn is_include_point2(&self, p: &[T; 2]) -> bool;
    fn corner_points(&self) -> [[T; 2]; 4];
    fn nearest_point2(&self, p: &[T; 2]) -> [T; 2];
    fn is_intersect_aabb2(&self, aabb: &[T; 4]) -> bool;
    fn is_intersect_obb2(&self, obb: &[T; 6]) -> bool;
}

impl OBB2Trait<f32> for [f32; 6] {
    fn is_include_point2(&self, p: &[f32; 2]) -> bool {
        is_include_point2(self, p)
    }
    fn corner_points(&self) -> [[f32; 2]; 4] {
        corner_points(self)
    }
    fn nearest_point2(&self, p: &[f32; 2]) -> [f32; 2] {
        nearest_point2(self, p)
    }
    fn is_intersect_aabb2(&self, aabb: &[f32; 4]) -> bool {
        is_intersect_aabb2(self, aabb)
    }
    fn is_intersect_obb2(&self, obb: &[f32; 6]) -> bool {
        is_intersect_obb2(self, obb)
    }
}

pub fn from_random<RAND>(reng: &mut RAND) -> [f32; 6]
where
    RAND: rand::Rng,
{
    let cntr = [2. * reng.gen::<f32>() - 1., 2. * reng.gen::<f32>() - 1.];
    let u = [2. * reng.gen::<f32>() - 1., 2. * reng.gen::<f32>() - 1.];
    let v = [2. * reng.gen::<f32>() - 1., 2. * reng.gen::<f32>() - 1.];
    let v = crate::vec2::orthogonalize(&u, &v);
    [cntr[0], cntr[1], u[0], u[1], v[0], v[1]]
}

fn is_include_point2(obb: &[f32; 6], p: &[f32; 2]) -> bool {
    let d = p.sub(obb[..2].try_into().unwrap());
    {
        let v = obb[2..4].try_into().unwrap();
        let vd = d.dot(v);
        let vv = v.dot(v);
        if vd < -vv || vd > vv {
            return false;
        }
    }
    {
        let v = obb[4..6].try_into().unwrap();
        let vd = d.dot(v);
        let vv = v.dot(v);
        if vd < -vv || vd > vv {
            return false;
        }
    }
    true
}

/// four corner points of obb2
/// sequentially connecting the outputs give a quadrilateral
pub fn corner_points(obb: &[f32; 6]) -> [[f32; 2]; 4] {
    [
        [obb[0] + obb[2] + obb[4], obb[1] + obb[3] + obb[5]],
        [obb[0] - obb[2] + obb[4], obb[1] - obb[3] + obb[5]],
        [obb[0] - obb[2] - obb[4], obb[1] - obb[3] - obb[5]],
        [obb[0] + obb[2] - obb[4], obb[1] + obb[3] - obb[5]],
    ]
}

pub fn nearest_point2(obb: &[f32; 6], p: &[f32; 2]) -> [f32; 2] {
    if is_include_point2(obb, p) {
        return *p;
    }
    let cp = corner_points(obb);
    let ps = [
        crate::edge2::nearest_point2(&cp[0], &cp[1], p).1,
        crate::edge2::nearest_point2(&cp[1], &cp[2], p).1,
        crate::edge2::nearest_point2(&cp[2], &cp[3], p).1,
        crate::edge2::nearest_point2(&cp[3], &cp[0], p).1,
    ];
    let ls = [
        crate::edge2::length(&ps[0], p),
        crate::edge2::length(&ps[1], p),
        crate::edge2::length(&ps[2], p),
        crate::edge2::length(&ps[3], p),
    ];
    let (i, _l) = ls
        .iter()
        .enumerate()
        .min_by(|&a, &b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    ps[i]
}

#[test]
fn test_nearest_point2() {
    let obb = [0., 0., 1., 1., -0.5, 0.5];
    let p0 = nearest_point2(&obb, &[1.5, 1.5]);
    assert!(crate::edge2::length(&p0, &[1., 1.]) < 1.0e-8);
}

/// check if an OBB is intersect with an AABB. For obb,
/// the first two element is one of its point, and the other are two edge vectors
pub fn is_intersect_aabb2(obb: &[f32; 6], aabb: &[f32; 4]) -> bool {
    // choose 4 separating axes for obb and aabb
    let axes = [[1.0, 0.0], [0.0, 1.0], [-obb[3], obb[2]], [-obb[5], obb[4]]];
    let obb_points = corner_points(obb);
    let aabb_points = [
        [aabb[0], aabb[1]],
        [aabb[2], aabb[1]],
        [aabb[0], aabb[3]],
        [aabb[2], aabb[3]],
    ];
    for axis in axes.iter() {
        let mut obb_min = f32::INFINITY;
        let mut obb_max = f32::NEG_INFINITY;
        let mut aabb_min = f32::INFINITY;
        let mut aabb_max = f32::NEG_INFINITY;
        for i in 0..4 {
            let obb_proj = crate::vec2::dot(&obb_points[i], axis);
            obb_min = obb_min.min(obb_proj);
            obb_max = obb_max.max(obb_proj);
            let aabb_proj = crate::vec2::dot(&aabb_points[i], axis);
            aabb_min = aabb_min.min(aabb_proj);
            aabb_max = aabb_max.max(aabb_proj);
        }
        if obb_max < aabb_min || obb_min > aabb_max {
            return false;
        }
    }
    true
}

#[test]
fn test_is_intersect_aabb2() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _iter in 0..100 {
        let obb = from_random(&mut reng);
        let aabb = crate::aabb2::from_two_points(
            &[2. * reng.gen::<f32>() - 1., 2. * reng.gen::<f32>() - 1.],
            &[2. * reng.gen::<f32>() - 1., 2. * reng.gen::<f32>() - 1.],
            0.,
        );
        assert!(aabb[0] < aabb[2]);
        assert!(aabb[1] < aabb[3]);
        let a0 = crate::aabb2::center(&aabb);
        let a1 = nearest_point2(&obb, &a0);
        let a2 = crate::aabb2::nearest_point2(&aabb, &a1);
        let a3 = nearest_point2(&obb, &a2);
        let a4 = crate::aabb2::nearest_point2(&aabb, &a3);
        let len23 = crate::edge2::length(&a2, &a3);
        let len34 = crate::edge2::length(&a3, &a4);
        if len34 > 0. && len34 < len23 * 0.99999 {
            continue;
        } // still converging
        let res0 = is_intersect_aabb2(&obb, &aabb);
        let res1 = len34 < 0.0001;
        assert_eq!(res0, res1);
    }
}

pub fn is_intersect_obb2(obb1: &[f32; 6], obb2: &[f32; 6]) -> bool {
    let axes = [
        [-obb1[3], obb1[2]],
        [-obb1[5], obb1[4]],
        [-obb2[3], obb2[2]],
        [-obb2[5], obb2[4]],
    ];
    let obb1_points = corner_points(obb1);
    let obb2_points = corner_points(obb2);
    for axis in axes.iter() {
        let mut obb1_min = f32::INFINITY;
        let mut obb1_max = f32::NEG_INFINITY;
        let mut obb2_min = f32::INFINITY;
        let mut obb2_max = f32::NEG_INFINITY;
        for i in 0..4 {
            let obb1_proj = obb1_points[i][0] * axis[0] + obb1_points[i][1] * axis[1];
            obb1_min = obb1_min.min(obb1_proj);
            obb1_max = obb1_max.max(obb1_proj);
            let obb2_proj = obb2_points[i][0] * axis[0] + obb2_points[i][1] * axis[1];
            obb2_min = obb2_min.min(obb2_proj);
            obb2_max = obb2_max.max(obb2_proj);
        }
        if obb1_max < obb2_min || obb1_min > obb2_max {
            return false;
        }
    }
    true
}

#[test]
fn test_is_intersect_obb1() {
    {
        let obb1 = [0., 0., 1.0, 2.0, 2.0, 1.0];
        let obb2 = [1., 1., 1.0, 2.0, 2.0, 1.0];
        assert!(is_intersect_obb2(&obb1, &obb2));
    }
    {
        let obb1 = [0., 0., 1.0, 0.0, 2.0, 1.0];
        let obb2 = [1.1, 0.0, 1.0, 0.0, 2.0, 1.0];
        assert!(is_intersect_obb2(&obb1, &obb2));
    }
}
