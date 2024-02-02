//! methods for 3x3 matrix

use num_traits::AsPrimitive;

pub fn minimum_rotation_matrix<T>(
    v0: nalgebra::Vector3::<T>,
    v1: nalgebra::Vector3::<T>) -> nalgebra::Matrix3::<T>
where T: nalgebra::RealField + 'static + Copy,
    f64: num_traits::AsPrimitive<T>
{
    use crate::vec3::frame_from_z_vector;

    let ep = v0.normalize();
    let eq = v1.normalize();
    let n = ep.cross(&eq);
    let st2 = n.dot(&n);
    let ct = ep.dot(&eq);
    let half = 0.5_f64.as_();

    if st2 < 1.0e-8_f64.as_() { // very small angle or n is zero
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
                eqx.dot(&epx), eqy.dot(&epx), eq.dot(&epx),
                eqx.dot(&epy), eqy.dot(&epy), eq.dot(&epy),
                eqx.dot(&ep), eqy.dot(&ep), eq.dot(&ep),
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
    is_increasing: bool) -> (nalgebra::Vector3<T>, nalgebra::Matrix3<T>)
where T: nalgebra::RealField + Copy
{
    let sgn = if is_increasing { T::one() } else {-T::one()};
    // let cov = evec*nalgebra::Matrix3::<T>::from_diagonal(eval)*evec.transpose();
    let mut prmt: Vec<usize> = vec!(0,1,2);
    prmt.sort_by(|&idx, &jdx| (sgn*eval[idx]).partial_cmp(&(sgn*eval[jdx])).unwrap() );
    // dbg!(&prmt,eval);
    let eval1 = nalgebra::Vector3::<T>::new(eval[prmt[0]], eval[prmt[1]], eval[prmt[2]]);
    let evec1 = nalgebra::Matrix3::<T>::from_columns(&[
        evec.column(prmt[0]), evec.column(prmt[1]), evec.column(prmt[2])]);
    //let cov1 = evec1*nalgebra::Matrix3::<T>::from_diagonal(&eval1)*evec1.transpose();
    (eval1, evec1)
}


/// when SVD of 3x3 matrix a is U*S*V^T, compute U*V^T
pub fn rotational_component<T>(
    a: &nalgebra::Matrix3::<T>) -> nalgebra::Matrix3::<T>
where T: nalgebra::RealField + Copy
{
    let svd = nalgebra::linalg::SVD::<T, nalgebra::U3, nalgebra::U3>::new(*a, true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let u_vt = u * v_t;
    let u_vt = if u_vt.determinant() > T::zero() { u_vt } else {
        let mut v_t = v_t;
        v_t.row_mut(0).scale_mut(-T::one());
        u * v_t
    };
    u_vt
}