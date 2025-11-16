use num_traits::AsPrimitive;
use num_traits::FloatConst;

pub fn to_uvec3_equal_area<Real>(p: &[Real; 2]) -> [Real; 3]
where
    Real: num_traits::Float + 'static + AsPrimitive<f64> + FloatConst,
    f64: AsPrimitive<Real>,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    debug_assert!(p[0] >= zero && p[0] <= one && p[1] >= zero && p[1] <= one);

    // Transform p to [-1, 1]^2 and compute absolute values
    let u = two * p[0] - one;
    let v = two * p[1] - one;
    let up = u.abs();
    let vp = v.abs();

    // Compute radius r as signed distance from diagonal
    let signed_distance = one - (up + vp);
    let d = signed_distance.abs();
    let r = one - d;

    // Compute angle φ for square to sphere mapping
    let phi = if r == zero { one } else { (vp - up) / r + one } * Real::FRAC_PI_4();

    let z = (one - r * r).copysign(signed_distance);

    // SafeSqrt(2 - r^2)
    let tmp = two - r * r;
    let safe_sqrt = if tmp > zero { tmp.sqrt() } else { zero };

    // Compute cosφ and sinφ for original quadrant and return vector
    let cos_phi = phi.cos().copysign(u);
    let sin_phi = phi.sin().copysign(v);

    [cos_phi * r * safe_sqrt, sin_phi * r * safe_sqrt, z]
}

#[test]
fn test_unit2_uvec3_unit2() {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0);
    for _itr in 0..1000 {
        let p0 = [rng.random::<f64>(), rng.random::<f64>()];
        let q = to_uvec3_equal_area(&p0);
        let p1 = crate::uvec3::map_to_unit2_equal_area(&q);
        use crate::vec2::Vec2;
        let diff = p0.sub(&p1).norm();
        assert!(diff < 1.0e-13, "{}", diff);
    }
}
