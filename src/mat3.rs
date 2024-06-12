//! methods for 3x3 matrix

use num_traits::AsPrimitive;

pub fn inverse_array_of_array<T>(gd: &[[T; 3]; 3]) -> [[T; 3]; 3]
where
    T: num_traits::Float,
{
    // contravariant basis vectors
    let mut gu = [[T::zero(); 3]; 3];
    crate::vec3::cross_mut_(&mut gu[0], &gd[1], &gd[2]);
    let invtmp1 = T::one() / crate::vec3::dot_(&gu[0], &gd[0]);
    gu[0][0] = gu[0][0] * invtmp1;
    gu[0][1] = gu[0][1] * invtmp1;
    gu[0][2] = gu[0][2] * invtmp1;
    //
    crate::vec3::cross_mut_(&mut gu[1], &gd[2], &gd[0]);
    let invtmp2 = T::one() / crate::vec3::dot_(&gu[1], &gd[1]);
    gu[1][0] = gu[1][0] * invtmp2;
    gu[1][1] = gu[1][1] * invtmp2;
    gu[1][2] = gu[1][2] * invtmp2;
    //
    crate::vec3::cross_mut_(&mut gu[2], &gd[0], &gd[1]);
    let invtmp3 = T::one() / crate::vec3::dot_(&gu[2], &gd[2]);
    gu[2][0] = gu[2][0] * invtmp3;
    gu[2][1] = gu[2][1] * invtmp3;
    gu[2][2] = gu[2][2] * invtmp3;
    gu
}

pub fn inverse<T>(b: &[T; 9]) -> [T; 9]
where
    T: num_traits::Float,
{
    let det = b[0] * b[4] * b[8] + b[3] * b[7] * b[2] + b[6] * b[1] * b[5]
        - b[0] * b[7] * b[5]
        - b[6] * b[4] * b[2]
        - b[3] * b[1] * b[8];
    let inv_det = T::one() / det;
    [
        inv_det * (b[4] * b[8] - b[5] * b[7]),
        inv_det * (b[2] * b[7] - b[1] * b[8]),
        inv_det * (b[1] * b[5] - b[2] * b[4]),
        inv_det * (b[5] * b[6] - b[3] * b[8]),
        inv_det * (b[0] * b[8] - b[2] * b[6]),
        inv_det * (b[2] * b[3] - b[0] * b[5]),
        inv_det * (b[3] * b[7] - b[4] * b[6]),
        inv_det * (b[1] * b[6] - b[0] * b[7]),
        inv_det * (b[0] * b[4] - b[1] * b[3]),
    ]
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

// above: no dependency
// ------------------------------
// below: dependency with nalgebra

pub fn minimum_rotation_matrix<T>(
    v0: nalgebra::Vector3<T>,
    v1: nalgebra::Vector3<T>,
) -> nalgebra::Matrix3<T>
where
    T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>,
{
    use crate::vec3::frame_from_z_vector;

    let ep = v0.normalize();
    let eq = v1.normalize();
    let n = ep.cross(&eq);
    let st2 = n.dot(&n);
    let ct = ep.dot(&eq);
    let half = 0.5_f64.as_();

    if st2 < 1.0e-8_f64.as_() {
        // very small angle or n is zero
        // inifinitesimal rotation
        if ct > 0.99_f64.as_() {
            return nalgebra::Matrix3::<T>::new(
                T::one() + half * (n.x * n.x - st2),
                -n.z + half * (n.x * n.y),
                n.y + half * (n.x * n.z),
                n.z + half * (n.y * n.x),
                T::one() + half * (n.y * n.y - st2),
                -n.x + half * (n.y * n.z),
                -n.y + half * (n.z * n.x),
                n.x + half * (n.z * n.y),
                T::one() + half * (n.z * n.z - st2),
            );
        } else {
            let (epx, epy) = frame_from_z_vector(ep);
            let eqx = epx - eq.scale(eq.dot(&epx)); // vector orthogonal to eq
            let eqy = eq.cross(&eqx);
            return nalgebra::Matrix3::<T>::new(
                eqx.dot(&epx),
                eqy.dot(&epx),
                eq.dot(&epx),
                eqx.dot(&epy),
                eqy.dot(&epy),
                eq.dot(&epy),
                eqx.dot(&ep),
                eqy.dot(&ep),
                eq.dot(&ep),
            );
        }
    }
    let st = st2.sqrt();
    let n = n.normalize();
    // Rodoriguez's rotation formula
    nalgebra::Matrix3::<T>::new(
        ct + (T::one() - ct) * n.x * n.x,
        -n.z * st + (T::one() - ct) * n.x * n.y,
        n.y * st + (T::one() - ct) * n.x * n.z,
        n.z * st + (T::one() - ct) * n.y * n.x,
        ct + (T::one() - ct) * n.y * n.y,
        -n.x * st + (T::one() - ct) * n.y * n.z,
        -n.y * st + (T::one() - ct) * n.z * n.x,
        n.x * st + (T::one() - ct) * n.z * n.y,
        ct + (T::one() - ct) * n.z * n.z,
    )
}

/// sort eigen value and reorder the eigen vector
/// TODO: write some test
pub fn sort_eigen<T>(
    eval: &nalgebra::Vector3<T>,
    evec: &nalgebra::Matrix3<T>,
    is_increasing: bool,
) -> (nalgebra::Vector3<T>, nalgebra::Matrix3<T>)
where
    T: nalgebra::RealField + Copy,
{
    let sgn = if is_increasing { T::one() } else { -T::one() };
    // let cov = evec*nalgebra::Matrix3::<T>::from_diagonal(eval)*evec.transpose();
    let mut prmt: Vec<usize> = vec![0, 1, 2];
    prmt.sort_by(|&idx, &jdx| (sgn * eval[idx]).partial_cmp(&(sgn * eval[jdx])).unwrap());
    // dbg!(&prmt,eval);
    let eval1 = nalgebra::Vector3::<T>::new(eval[prmt[0]], eval[prmt[1]], eval[prmt[2]]);
    let evec1 = nalgebra::Matrix3::<T>::from_columns(&[
        evec.column(prmt[0]),
        evec.column(prmt[1]),
        evec.column(prmt[2]),
    ]);
    //let cov1 = evec1*nalgebra::Matrix3::<T>::from_diagonal(&eval1)*evec1.transpose();
    (eval1, evec1)
}

/// when SVD of 3x3 matrix a is U*S*V^T, compute U*V^T
pub fn rotational_component<T>(a: &nalgebra::Matrix3<T>) -> nalgebra::Matrix3<T>
where
    T: nalgebra::RealField + Copy,
{
    let svd = nalgebra::linalg::SVD::<T, nalgebra::U3, nalgebra::U3>::new(*a, true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let u_vt = u * v_t;
    let u_vt = if u_vt.determinant() > T::zero() {
        u_vt
    } else {
        let mut v_t = v_t;
        v_t.row_mut(0).scale_mut(-T::one());
        u * v_t
    };
    u_vt
}

pub fn skew<T>(v: &nalgebra::Vector3<T>) -> nalgebra::Matrix3<T>
where
    T: nalgebra::RealField + Copy,
{
    nalgebra::Matrix3::new(
        T::zero(),
        -v[2],
        v[1],
        v[2],
        T::zero(),
        -v[0],
        -v[1],
        v[0],
        T::zero(),
    )
}

#[test]
fn test_skew() {
    let v0 = nalgebra::Vector3::<f64>::new(1.1, 3.1, 2.5);
    let v1 = nalgebra::Vector3::<f64>::new(2.1, 0.1, 4.5);
    let c0 = v0.cross(&v1);
    let c1 = skew(&v0) * v1;
    assert!((c0 - c1).norm() < 1.0e-10);
}

/*
pub fn to_na_from_row_major_slice<T>(s: &[T;9]) -> nalgebra::Matrix3::<T>
    where T: Copy + nalgebra::RealField
{
    nalgebra::Matrix3::<T>::new(
        s[0], s[3], s[6],
        s[1], s[4], s[7],
        s[2], s[5], s[8])
}
 */
