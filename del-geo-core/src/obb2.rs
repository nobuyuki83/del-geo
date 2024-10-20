//! 2D Oriented Bounding Box
//! data structure `&[Real;6]`
//! first 2 Reals are for center
//! next 2 Reals are for half of major axis direction
//! next 2 Reals are for half of minar axis direction

use core::f32;

use crate::aabb;


/// check if an OBB is intersect with an AABB
pub fn is_intersect_aabb_2d(obb: &[f32;6], aabb: &[f32;4]) -> bool {
    // choose 4 separating axes for obb and aabb
    let axes = [[1.0, 0.0], [0.0, 1.0], [obb[2], obb[3]], [obb[4], obb[5]]];
    let obb_points = [
        [obb[0]                  , obb[1]                  ],
        [obb[0] + obb[2]         , obb[1] + obb[3]         ],
        [obb[0] + obb[4]         , obb[1] + obb[5]         ],
        [obb[0] + obb[2] + obb[4], obb[1] + obb[3] + obb[5]],
    ];
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
            let obb_proj = obb_points[i][0] * axis[0] + obb_points[i][1] * axis[1];
            obb_min = obb_min.min(obb_proj);
            obb_max = obb_max.max(obb_proj);
            let aabb_proj = aabb_points[i][0] * axis[0] + aabb_points[i][1] * axis[1];
            aabb_min = aabb_min.min(aabb_proj);
            aabb_max = aabb_max.max(aabb_proj);
        }
        if obb_max < aabb_min || obb_min > aabb_max {
            return false;
        }
    }
    true
}

pub fn is_intersect_obb_2d(obb1: &[f32;6], obb2: &[f32;6]) -> bool {
    let axes = [[obb1[2], obb1[3]], [obb1[4], obb1[5]], [obb2[2], obb2[3]], [obb2[4], obb2[5]]];
    let obb1_points = [
        [obb1[0]                    , obb1[1]                  ],
        [obb1[0] + obb1[2]          , obb1[1] + obb1[3]         ],
        [obb1[0] + obb1[4]          , obb1[1] + obb1[5]         ],
        [obb1[0] + obb1[2] + obb1[4], obb1[1] + obb1[3] + obb1[5]],
    ];
    let obb2_points = [
        [obb2[0]                    , obb2[1]                  ],
        [obb2[0] + obb2[2]          , obb2[1] + obb2[3]         ],
        [obb2[0] + obb2[4]          , obb2[1] + obb2[5]         ],
        [obb2[0] + obb2[2] + obb2[4], obb2[1] + obb2[3] + obb2[5]],
    ];
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
fn test_is_intersect_aabb1() {
    let obb = [0., 0., 1.0, 0.0, 0.0, 1.0];
    let aabb = [0.5, 0.5, 1.5, 1.5];
    assert!(is_intersect_aabb_2d(&obb, &aabb));
}

#[test]
fn test_is_intersect_aabb2() {
    let obb = [0., 0., 1.0, 1.0, -1.0, 1.0];
    let aabb = [-0.5, -0.5, -1.5, -1.5];
    assert!(!is_intersect_aabb_2d(&obb, &aabb));
}

#[test]
fn test_is_intersect_obb() {
    let obb1 = [0., 0., 1.0, 2.0, 2.0, 1.0];
    let obb2 = [1., 1., 1.0, 2.0, 2.0, 1.0];
    assert!(is_intersect_obb_2d(&obb1, &obb2));
}