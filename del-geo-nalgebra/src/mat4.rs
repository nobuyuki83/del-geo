pub fn jacobian_transform<Real>(
    t: &nalgebra::Matrix4::<Real>,
    p: &nalgebra::Vector3::<Real>) -> nalgebra::Matrix3<Real>
    where Real: nalgebra::RealField + Copy
{
    let a = nalgebra::Matrix3::<Real>::new(
        t.m11, t.m12, t.m13,
        t.m21, t.m22, t.m23,
        t.m31, t.m32, t.m33);
    let b = nalgebra::Vector3::new(t.m14, t.m24, t.m34);
    let d = t.m44;
    let c = nalgebra::Vector3::new(t.m41, t.m42, t.m43);
    let cpd = c.dot(p) + d;
    let apb: nalgebra::Vector3<Real> = a * p + b;
    a / cpd - apb * c.transpose() / (cpd * cpd)
}