/// methods for 3x3 matrix

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