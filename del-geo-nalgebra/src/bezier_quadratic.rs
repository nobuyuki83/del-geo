//! functions for quadratic Bezier curve

pub fn eval<Real, const N: usize>(
    p0: &[Real; N],
    p1: &[Real; N],
    p2: &[Real; N],
    t0: Real,
) -> [Real; N]
where
    Real: num_traits::Float + Copy + std::iter::Sum,
{
    let one = Real::one();
    let two = one + one;
    let t1 = one - t0;
    use del_geo_core::vecn::VecN;
    del_geo_core::vecn::add_three_vectors(
        &p0.scale(t1 * t1),
        &p1.scale(two * t0 * t1),
        &p2.scale(t0 * t0),
    )
}
