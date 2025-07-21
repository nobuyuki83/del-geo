pub fn mult_vec<Real, const N: usize, const NN: usize>(a: &[Real; NN], b: &[Real; N]) -> [Real; N]
where
    Real: num_traits::Float,
{
    assert_eq!(NN, N * N);
    let mut c = [Real::zero(); N];
    for i in 0..N {
        for j in 0..N {
            c[i] = c[i] + a[i + N * j] * b[j];
        }
    }
    c
}

pub fn add_in_place<Real, const NN: usize>(a: &mut [Real; NN], b: &[Real; NN])
where
    Real: num_traits::Float,
{
    a.iter_mut()
        .zip(b.iter())
        .for_each(|(va, vb)| *va = *va + *vb);
}

pub fn mult_mat_col_major<Real, const N: usize, const NN: usize>(
    a: &[Real; NN],
    b: &[Real; NN],
) -> [Real; NN]
where
    Real: num_traits::Float,
{
    assert_eq!(NN, N * N);
    let mut c = [Real::zero(); NN];
    for i in 0..N {
        for j in 0..N {
            for k in 0..N {
                c[i + N * j] = c[i + N * j] + a[i + N * k] * b[k + N * j];
            }
        }
    }
    c
}

pub fn try_inverse<Real, const N: usize, const NN: usize>(b: &[Real; NN]) -> Option<[Real; NN]>
where
    Real: num_traits::Float,
{
    crate::matn_row_major::try_inverse::<Real, N, NN>(b)
}

pub fn sub_in_place<Real, const NN: usize>(a: &mut [Real; NN], b: &[Real; NN])
where
    Real: num_traits::Float,
{
    a.iter_mut()
        .zip(b.iter())
        .for_each(|(va, &vb)| *va = (*va) - vb);
}

pub fn transpose<Real, const N: usize, const NN: usize>(a: &[Real; NN]) -> [Real; NN]
where
    Real: num_traits::Float,
{
    let mut c = [Real::zero(); NN];
    for i in 0..N {
        for j in 0..N {
            c[i + N * j] = a[j + N * i];
        }
    }
    c
}
