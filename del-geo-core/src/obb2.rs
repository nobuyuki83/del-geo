//! 2D Oriented Bounding Box
//! data structure `&[Real;6]`
//! first 2 Reals are for center
//! next 2 Reals are for half of major axis direction
//! next 2 Reals are for half of minar axis direction

pub fn is_intersect_aabb(obb: &[f32;6], aabb: &[f32;4]) -> bool {
    // hint use SAT (Separating Axis Theorem)
    todo!();
}


#[test]
fn test_is_intersect_aabb() {
    let obb = [0., 0., 1.0, 0.0, 0.0, 1.0];
    let aabb = [0.5, 0.5, 1.5, 1.5];
    // assert!(is_intersect_aabb(&obb, &aabb));
}

