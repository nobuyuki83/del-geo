//! methods for 3x3 matrix

pub trait Mat3ArrayOfArray<Real>
where
    Self: Sized,
{
    fn det_inv(&self) -> (Real, Self);
    fn inverse(&self) -> Self;
    fn matmul(&self, b: &Self) -> Self;
}

impl<Real> Mat3ArrayOfArray<Real> for [[Real; 3]; 3]
where
    Real: num_traits::Float,
{
    fn det_inv(&self) -> (Real, Self) {
        det_inv(self)
    }
    fn inverse(&self) -> Self {
        inverse(self)
    }
    fn matmul(&self, b: &Self) -> Self {
        matmul(self, b)
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
    a.det_inv().1
}

pub fn matmul<T>(a: &[[T; 3]; 3], b: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    let mut result = [[T::zero(); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            #[allow(clippy::needless_range_loop)]
            for k in 0..3 {
                result[i][j] = result[i][j] + a[i][k] * b[k][j];
            }
        }
    }
    result
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

pub fn from_identity<T>() -> [[T;3];3]
where T: num_traits::Float
{
    let zero = T::zero();
    let one = T::one();
    [
        [one, zero, zero],
        [zero, one, zero],
        [zero, zero, one],
    ]
}