//! methods for 3d triangle
use std::ops::MulAssign;

use crate::vec3::Vec3;

/// clamp barycentric coordinates inside a triangle
pub fn clamp<T>(r0: T, r1: T, r2: T) -> (T, T, T)
where
    T: num_traits::Float,
{
    // vertex projection
    if r0 <= T::zero() && r1 <= T::zero() {
        return (T::zero(), T::zero(), T::one());
    }
    if r1 <= T::zero() && r2 <= T::zero() {
        return (T::one(), T::zero(), T::zero());
    }
    if r2 <= T::zero() && r0 <= T::zero() {
        return (T::zero(), T::one(), T::zero());
    }
    // edge projection
    if r0 <= T::zero() {
        return (T::zero(), r1 / (r1 + r2), r2 / (r1 + r2));
    }
    if r1 <= T::zero() {
        return (r0 / (r0 + r2), T::zero(), r2 / (r0 + r2));
    }
    if r2 <= T::zero() {
        return (r0 / (r0 + r1), r1 / (r0 + r1), T::zero());
    }
    (r0, r1, r2)
}

// ----------------------

/// normal vector of a 3D triangle (coordinates given by stack-allocated arrays)
pub fn normal<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3]) -> [T; 3]
where
    T: std::ops::Sub<Output = T> + std::ops::Mul<Output = T> + std::ops::Sub + Copy,
{
    [
        (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
        (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
        (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]),
    ]
}

/// area of a 3D triangle (coordinates given by stack-allocated arrays)
pub fn area<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3]) -> T
where
    T: num_traits::Float + MulAssign,
{
    let na = normal(v1, v2, v3);
    let half = T::one() / (T::one() + T::one());
    na.squared_norm().sqrt() * half
}

pub fn unit_normal_area<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> ([T; 3], T)
where
    T: num_traits::Float + MulAssign,
{
    let n = normal(p0, p1, p2);
    let half = T::one() / (T::one() + T::one());
    let a = n.norm() * half;
    let invlen: T = half / a;
    ([n[0] * invlen, n[1] * invlen, n[2] * invlen], a)
}

/// compute cotangents of the three angles of a triangle
pub fn cot<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float + MulAssign,
{
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let v1 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let v2 = [p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]];
    let half = T::one() / (T::one() + T::one());
    let onefourth = T::one() / (T::one() + T::one() + T::one() + T::one());
    let area = {
        let na = [
            v1[1] * v2[2] - v2[1] * v1[2],
            v1[2] * v2[0] - v2[2] * v1[0],
            v1[0] * v2[1] - v2[0] * v1[1],
        ];
        na.squared_norm().sqrt() * half
    };
    let tmp: T = onefourth / area;
    let l0 = v0.squared_norm();
    let l1 = v1.squared_norm();
    let l2 = v2.squared_norm();
    [
        (l1 + l2 - l0) * tmp, // cot0 = cos0/sin0 = {(l1*l1+l2*l2-l0*l0)/(2*l1*l2)} / {2*area/(l1*l2)}
        (l2 + l0 - l1) * tmp, // cot1 = cos1/sin1 = {(l2*l2+l0*l0-l1*l1)/(2*l2*l0)} / {2*area/(l2*l0)}
        (l0 + l1 - l2) * tmp, // cot2 = cos2/sin2 = {(l0*l0+l1*l1-l2*l2)/(2*l0*l1)} / {2*area/(l0*l1)}
    ]
}

pub fn emat_cotangent_laplacian<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [[[T; 1]; 3]; 3]
where
    T: num_traits::Float + MulAssign,
{
    let cots = cot(p0, p1, p2);
    let half = T::one() / (T::one() + T::one());
    [
        [
            [(cots[1] + cots[2]) * half],
            [-cots[2] * half],
            [-cots[1] * half],
        ],
        [
            [-cots[2] * half],
            [(cots[2] + cots[0]) * half],
            [-cots[0] * half],
        ],
        [
            [-cots[1] * half],
            [-cots[0] * half],
            [(cots[0] + cots[1]) * half],
        ],
    ]
}

pub fn emat_graph_laplacian<T>(l: T) -> [[[T; 1]; 3]; 3]
where
    T: num_traits::Float,
{
    let vo = -T::one() * l;
    let vd = (T::one() + T::one()) * l;
    [[[vd], [vo], [vo]], [[vo], [vd], [vo]], [[vo], [vo], [vd]]]
}

/// Möller–Trumbore ray-triangle intersection algorithm
///
/// <https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm>
/// <https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection>
pub fn intersection_against_line<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    ray_org: &[T; 3],
    ray_dir: &[T; 3],
) -> Option<T>
where
    T: num_traits::Float + MulAssign,
{
    let eps: T = T::epsilon();
    let edge1 = p1.sub(p0);
    let edge2 = p2.sub(p0);
    let pvec = ray_dir.cross(&edge2);
    let det = edge1.dot(&pvec);
    if det > -eps && det < eps {
        return None;
    }
    let invdet = T::one() / det;
    let tvec = ray_org.sub(p0);
    let u = invdet * tvec.dot(&pvec);
    if u < T::zero() || u > T::one() {
        return None;
    }
    let qvec = tvec.cross(&edge1);
    let v = invdet * ray_dir.dot(&qvec);
    if v < T::zero() || u + v > T::one() {
        return None;
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    let t = invdet * edge2.dot(&qvec);
    Some(t)
}

pub struct RayTriangleIntersectionData<Real> {
    inv_det: Real,
    n_dot_c: Real,
    n_dot_dir: Real,
    r_dot_e2: Real,
    r_dot_e1: Real,
    c: [Real; 3],
    n: [Real; 3],
    r: [Real; 3],
    e1: [Real; 3],
    e2: [Real; 3],
    pub dir: [Real; 3],
}

/// ray triangle intersection.
/// * `dir` - any nonzero vector (not necessary to be a unit vector)
/// * `t` - ratio of `dir` vector from
/// * `u` - barycentric coordinate
/// * `v` - barycentric coordinate
///
/// `org + t * dir = (1 - u - v) * p0 + u * p1 + v * p2`
pub fn ray_triangle_intersection<T>(
    org: &[T; 3],
    dir: &[T; 3],
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
) -> Option<(T, T, T, RayTriangleIntersectionData<T>)>
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let e1 = p0.sub(p1);
    let e2 = p2.sub(p0);
    let n = e1.cross(&e2);
    let n_dot_dir = n.dot(dir);
    if n_dot_dir.is_zero() {
        return None;
    }
    let inv_det = T::one() / n_dot_dir;
    let c = p0.sub(org);
    let n_dot_c = n.dot(&c);
    let t_ = inv_det * n_dot_c;
    if t_ < T::zero() {
        return None;
    }
    let r = dir.cross(&c);
    //
    let r_dot_e2 = r.dot(&e2);
    let u_ = inv_det * r_dot_e2;
    if u_ < T::zero() {
        return None;
    }
    //
    let r_dot_e1 = r.dot(&e1);
    let v_ = inv_det * r_dot_e1;
    if v_ < T::zero() {
        return None;
    }

    if u_ + v_ >= T::one() {
        return None;
    }
    Some((
        t_,
        u_,
        v_,
        RayTriangleIntersectionData {
            inv_det,
            n_dot_c,
            n_dot_dir,
            r_dot_e2,
            r_dot_e1,
            c,
            n,
            r,
            e1,
            e2,
            dir: *dir,
        },
    ))
}

#[test]
fn test_ray_triangle_intersection() {
    let p0 = [-13f32, -5., 8.];
    let p1 = [14., -5., 8.];
    let p2 = [1., 3., -3.];
    let origin = [8., 11., 10.];
    let dir = [-7., -11., -8.];
    let t1 = intersection_against_line(&p0, &p1, &p2, &origin, &dir).unwrap();
    let (t0, _u0, _v0, _data) = ray_triangle_intersection(&origin, &dir, &p0, &p1, &p2).unwrap();
    assert!((t0 - t1).abs() < 1.0e-8);
}

pub fn dldw_ray_triangle_intersection_<Real>(
    d_t: Real,
    d_u: Real,
    d_v: Real,
    data: &RayTriangleIntersectionData<Real>,
) -> ([Real; 3], [Real; 3], [Real; 3])
where
    Real: num_traits::Float + Copy,
{
    let d_n_dot_c = d_t * data.inv_det;
    let d_r_dot_e2 = d_u * data.inv_det;
    let d_r_dot_e1 = d_v * data.inv_det;
    let d_inv_det = d_t * data.n_dot_c + d_u * data.r_dot_e2 + d_v * data.r_dot_e1;
    //
    let mut d_n = data.c.scale(d_n_dot_c);
    let mut d_c = data.n.scale(d_n_dot_c);
    let mut d_e2 = data.r.scale(d_r_dot_e2);
    let mut d_e1 = data.r.scale(d_r_dot_e1);
    let d_r = data.e2.scale(d_r_dot_e2).add(&data.e1.scale(d_r_dot_e1));
    //
    let d_n_dot_dir = -d_inv_det / data.n_dot_dir / data.n_dot_dir;
    d_n.add_in_place(&data.dir.scale(d_n_dot_dir));
    d_c.add_in_place(&d_r.cross(&data.dir));
    d_e2.add_in_place(&d_n.cross(&data.e1));
    d_e1.add_in_place(&data.e2.cross(&d_n));
    //
    let d_p0 = d_e1.sub(&d_e2).add(&d_c);
    let d_p1 = d_e1.scale(-Real::one());
    let d_p2 = d_e2;
    (d_p0, d_p1, d_p2)
}

#[test]
fn test_dw_ray_triangle_intersection() {
    type Real = f64;
    let p0: [[Real; 3]; 3] = [[-13., -5., 8.], [14., -5., 8.], [1., 3., -3.]];

    // let origin = nalgebra::Vector3::<Real>::new(9., 11., 12.);
    let origin = [8., 11., 10.];
    let dir = [1., 0., 2.].sub(&origin);

    let Some((t0, u0, v0, data)) = ray_triangle_intersection(&origin, &dir, &p0[0], &p0[1], &p0[2])
    else {
        panic!()
    };

    let ha = p0[0]
        .scale(1. - u0 - v0)
        .add(&p0[1].scale(u0))
        .add(&p0[2].scale(v0));
    assert!(crate::tet::volume(&p0[0], &p0[1], &p0[2], &ha).abs() < 1.0e-8);
    let hb = crate::vec3::axpy(t0, &dir, &origin);
    assert!(ha.sub(&hb).norm() < 1.0e-6);

    // d_t, d_u, d_u are back-propagated from the loss
    let d_t = 0.7;
    let d_u = 1.3;
    let d_v = 1.1;

    let dp = dldw_ray_triangle_intersection_::<Real>(d_t, d_u, d_v, &data);
    let dp = [dp.0, dp.1, dp.2];

    let eps = 1.0e-5;
    for i_node in 0..3 {
        for i_dim in 0..3 {
            let p1 = {
                let mut p1 = p0;
                p1[i_node][i_dim] += eps;
                p1
            };
            let Some((t1, u1, v1, _data)) =
                ray_triangle_intersection(&origin, &dir, &p1[0], &p1[1], &p1[2])
            else {
                panic!()
            };
            let dloss = ((t1 - t0) * d_t + (u1 - u0) * d_u + (v1 - v0) * d_v) / eps;
            let dloss_analytic = dp[i_node][i_dim];
            // dbg!(dloss, dloss_analytic);
            assert!(
                (dloss - dloss_analytic).abs() < 6.0e-7,
                "{} {} {}",
                dloss_analytic,
                dloss,
                dloss - dloss_analytic
            );
        }
    }
}

// ----------------------------

pub fn to_barycentric_coords<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3], q: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float + MulAssign,
{
    let a0 = area(q, p1, p2);
    let a1 = area(q, p2, p0);
    let a2 = area(q, p0, p1);
    let sum_inv = T::one() / (a0 + a1 + a2);
    [a0 * sum_inv, a1 * sum_inv, a2 * sum_inv]
}

pub fn position_from_barycentric_coords<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    bc: &[T; 3],
) -> [T; 3]
where
    T: num_traits::Float,
{
    [
        bc[0] * p0[0] + bc[1] * p1[0] + bc[2] * p2[0],
        bc[0] * p0[1] + bc[1] * p1[1] + bc[2] * p2[1],
        bc[0] * p0[2] + bc[1] * p1[2] + bc[2] * p2[2],
    ]
}

// -------------------------

pub struct Tri3<'a, Real> {
    pub p0: &'a [Real; 3],
    pub p1: &'a [Real; 3],
    pub p2: &'a [Real; 3],
}

#[allow(clippy::needless_lifetimes)]
impl<'a, Real> Tri3<'a, Real>
where
    Real: num_traits::Float + MulAssign,
{
    pub fn intersection_against_ray(
        &self,
        ray_org: &[Real; 3],
        ray_dir: &[Real; 3],
    ) -> Option<Real> {
        intersection_against_line(self.p0, self.p1, self.p2, ray_org, ray_dir)
            .filter(|&t| t >= Real::zero())
    }

    pub fn intersection_against_line(
        &self,
        line_org: &[Real; 3],
        line_dir: &[Real; 3],
    ) -> Option<Real> {
        intersection_against_line(self.p0, self.p1, self.p2, line_org, line_dir)
    }

    pub fn area(&self) -> Real {
        area(self.p0, self.p1, self.p2)
    }

    pub fn cog(&self) -> [Real; 3] {
        let one = Real::one();
        let one3rd = one / (one * one * one);
        [
            (self.p0[0] + self.p1[0] + self.p2[0]) * one3rd,
            (self.p0[1] + self.p1[1] + self.p2[1]) * one3rd,
            (self.p0[2] + self.p1[2] + self.p2[2]) * one3rd,
        ]
    }

    pub fn normal(&self) -> [Real; 3] {
        normal(self.p0, self.p1, self.p2)
    }

    pub fn unit_normal(&self) -> [Real; 3] {
        let n = normal(self.p0, self.p1, self.p2);
        n.normalize()
    }

    pub fn position_from_barycentric_coordinates(&self, r0: Real, r1: Real) -> [Real; 3] {
        let r2 = Real::one() - r0 - r1;
        [
            r0 * self.p0[0] + r1 * self.p1[0] + r2 * self.p2[0],
            r0 * self.p0[1] + r1 * self.p1[1] + r2 * self.p2[1],
            r0 * self.p0[2] + r1 * self.p1[2] + r2 * self.p2[2],
        ]
    }
}
