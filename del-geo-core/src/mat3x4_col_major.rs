pub fn from_mat4_col_major<Real>(a: &[Real; 16]) -> [Real; 12]
where
    Real: num_traits::Float,
{
    [
        a[0], a[1], a[2], a[4], a[5], a[6], a[8], a[9], a[10], a[12], a[13], a[14],
    ]
}

pub fn transform_affine<Real>(a: &[Real; 12], v: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    [
        a[0] * v[0] + a[3] * v[1] + a[6] * v[2] + a[9],
        a[1] * v[0] + a[4] * v[1] + a[7] * v[2] + a[10],
        a[2] * v[0] + a[5] * v[1] + a[8] * v[2] + a[11],
    ]
}
