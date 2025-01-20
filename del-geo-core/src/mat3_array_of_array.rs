//! methods for 3x3 matrix

pub trait Mat3ArrayOfArray<Real>
where
    Self: Sized,
{
    fn det_inv(self) -> (Real, Self);
    fn inverse(self) -> Self;
    fn matmul(self, b: &Self) -> Self;
}

impl<Real> Mat3ArrayOfArray<Real> for [[Real; 3]; 3]
where
    Real: num_traits::Float,
{
    fn det_inv(self) -> (Real, Self) {
        det_inv(&self)
    }
    fn inverse(self) -> Self {
        inverse(&self)
    }
    fn matmul(self, b: &Self) -> Self {
        matmul(&self, b)
    }
}

pub fn det_inv<T>(a: &[[T; 3]; 3]) -> (T, [[T; 3]; 3])
where
    T: num_traits::Float,
{
    let det =
        a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
            - a[0][0] * a[2][1] * a[1][2]
            - a[2][0] * a[1][1] * a[0][2]
            - a[1][0] * a[0][1] * a[2][2];
    let inv_det = T::one() / det;
    let ainv = [
        [
            inv_det * (a[1][1] * a[2][2] - a[1][2] * a[2][1]),
            inv_det * (a[0][2] * a[2][1] - a[0][1] * a[2][2]),
            inv_det * (a[0][1] * a[1][2] - a[0][2] * a[1][1]),
        ],
        [
            inv_det * (a[1][2] * a[2][0] - a[1][0] * a[2][2]),
            inv_det * (a[0][0] * a[2][2] - a[0][2] * a[2][0]),
            inv_det * (a[0][2] * a[1][0] - a[0][0] * a[1][2]),
        ],
        [
            inv_det * (a[1][0] * a[2][1] - a[1][1] * a[2][0]),
            inv_det * (a[0][1] * a[2][0] - a[0][0] * a[2][1]),
            inv_det * (a[0][0] * a[1][1] - a[0][1] * a[1][0]),
        ],
    ];
    (det, ainv)
}

pub fn inverse<T>(a: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    det_inv(a).1
}

pub fn matmul<T>(a: &[[T; 3]; 3], b: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1],
            a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1],
            a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2],
        ],
        [
            a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0],
            a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1],
            a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2],
        ],
    ]
}

#[test]
fn test_inverse_matmul() {
    let a = [[0f64, 2., 4.], [3., 5., 4.], [6., 7., 8.]];
    let ainv = a.inverse();
    let ainv_a = ainv.matmul(&a);
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                assert!((1.0 - ainv_a[i][j]).abs() < f64::EPSILON);
            } else {
                assert!(ainv_a[i][j].abs() < f64::EPSILON);
            }
        }
    }
}
