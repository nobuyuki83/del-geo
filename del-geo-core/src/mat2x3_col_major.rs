

pub fn mult_vec3<Real>(a: &[Real;6], b: &[Real;9]) -> [Real;6]
where Real: num_traits::Float
{
    [
        a[0] * b[0] + a[2] * b[1] + a[4] * b[2],
        a[1] * b[0] + a[3] * b[1] + a[5] * b[2],
        a[0] * b[3] + a[2] * b[4] + a[4] * b[5],
        a[1] * b[3] + a[3] * b[4] + a[5] * b[5],
        a[0] * b[6] + a[2] * b[7] + a[4] * b[8],
        a[1] * b[6] + a[3] * b[7] + a[5] * b[8],
    ]
}
