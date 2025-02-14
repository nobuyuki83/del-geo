pub fn from_mat4_col_major<Real>(a: &[Real; 16]) -> [Real; 12]
where
    Real: num_traits::Float,
{
    [
        a[0], a[1], a[2], a[4], a[5], a[6], a[8], a[9], a[10], a[12], a[13], a[14],
    ]
}
