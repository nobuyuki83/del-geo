
pub fn eval<Real, const N: usize>(
    p0: &nalgebra::SVector<Real, N>,
    p1: &nalgebra::SVector<Real, N>,
    p2: &nalgebra::SVector<Real, N>,
    p3: &nalgebra::SVector<Real, N>,
    t0: Real) -> nalgebra::SVector<Real, N>
    where Real: nalgebra::RealField + Copy
{
    let one = Real::one();
    let three = one + one + one;
    let t1 = one - t0;

    p0.scale( t1 * t1 * t1)
        + p1.scale( three * t0 * t1 * t1)
        + p2.scale( three * t0 * t0 * t1 )
        + p3.scale( t0 * t0 * t0)
}