pub fn try_inverse<Real, const N: usize, const NN: usize>(b: &[Real; NN]) -> Option<[Real; NN]>
where
    Real: num_traits::Float + std::ops::MulAssign + std::ops::SubAssign,
{
    let mut a = *b;
    for i in 0..N {
        if a[i * N + i].is_zero() {
            return None;
        }
        {
            let tmp1 = Real::one() / a[i * N + i];
            a[i * N + i] = Real::one();
            for k in 0..N {
                a[i * N + k] *= tmp1;
            }
        }
        for j in 0..N {
            if j == i {
                continue;
            }
            let tmp2 = a[j * N + i];
            a[j * N + i] = Real::zero();
            for k in 0..N {
                a[j * N + k] -= tmp2 * a[i * N + k];
            }
        }
    }
    Some(a)
}
