//! methods for 3x3 matrix

pub fn inverse_array_of_array<T>(gd: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    // contravariant basis vectors
    let mut gu = [[T::zero(); 3]; 3];
    crate::vec3::cross_mut_(&mut gu[0], &gd[1], &gd[2]);
    let invtmp1 = T::one() / crate::vec3::dot(&gu[0], &gd[0]);
    gu[0][0] = gu[0][0] * invtmp1;
    gu[0][1] = gu[0][1] * invtmp1;
    gu[0][2] = gu[0][2] * invtmp1;
    //
    crate::vec3::cross_mut_(&mut gu[1], &gd[2], &gd[0]);
    let invtmp2 = T::one() / crate::vec3::dot(&gu[1], &gd[1]);
    gu[1][0] = gu[1][0] * invtmp2;
    gu[1][1] = gu[1][1] * invtmp2;
    gu[1][2] = gu[1][2] * invtmp2;
    //
    crate::vec3::cross_mut_(&mut gu[2], &gd[0], &gd[1]);
    let invtmp3 = T::one() / crate::vec3::dot(&gu[2], &gd[2]);
    gu[2][0] = gu[2][0] * invtmp3;
    gu[2][1] = gu[2][1] * invtmp3;
    gu[2][2] = gu[2][2] * invtmp3;
    gu
}

pub fn try_inverse<T>(b: &[T; 9]) -> Option<[T; 9]>
where
    T: num_traits::Float,
{
    let det = b[0] * b[4] * b[8] + b[3] * b[7] * b[2] + b[6] * b[1] * b[5]
        - b[0] * b[7] * b[5]
        - b[6] * b[4] * b[2]
        - b[3] * b[1] * b[8];
    if det.is_zero() {
        return None;
    }
    let inv_det = T::one() / det;
    Some([
        inv_det * (b[4] * b[8] - b[5] * b[7]),
        inv_det * (b[2] * b[7] - b[1] * b[8]),
        inv_det * (b[1] * b[5] - b[2] * b[4]),
        inv_det * (b[5] * b[6] - b[3] * b[8]),
        inv_det * (b[0] * b[8] - b[2] * b[6]),
        inv_det * (b[2] * b[3] - b[0] * b[5]),
        inv_det * (b[3] * b[7] - b[4] * b[6]),
        inv_det * (b[1] * b[6] - b[0] * b[7]),
        inv_det * (b[0] * b[4] - b[1] * b[3]),
    ])
}

pub fn transform_homogeneous<Real>(transform: &[Real; 9], x: &[Real; 2]) -> Option<[Real; 2]>
where
    Real: num_traits::Float,
{
    let y2 = transform[2] * x[0] + transform[5] * x[1] + transform[8];
    if y2.is_zero() {
        return None;
    }
    //
    let y0 = transform[0] * x[0] + transform[3] * x[1] + transform[6];
    let y1 = transform[1] * x[0] + transform[4] * x[1] + transform[7];
    Some([y0 / y2, y1 / y2])
}

pub fn identity<T>() -> [T; 9]
where
    T: num_traits::Float,
{
    let zero = T::zero();
    let one = T::one();
    [one, zero, zero, zero, one, zero, zero, zero, one]
}
