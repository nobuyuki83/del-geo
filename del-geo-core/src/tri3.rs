//! methods for 3d triangle

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
    T: num_traits::Float,
{
    let na = normal(v1, v2, v3);
    let half = T::one() / (T::one() + T::one());
    crate::vec3::squared_norm(&na).sqrt() * half
}

pub fn unit_normal_area<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> ([T; 3], T)
where
    T: num_traits::Float,
{
    use crate::vec3;
    let n = normal(p0, p1, p2);
    let half = T::one() / (T::one() + T::one());
    let a = vec3::norm(&n) * half;
    let invlen: T = half / a;
    ([n[0] * invlen, n[1] * invlen, n[2] * invlen], a)
}

/// compute cotangents of the three angles of a triangle
pub fn cot<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    use crate::vec3;
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
        vec3::squared_norm(&na).sqrt() * half
    };
    let tmp: T = onefourth / area;
    let l0 = vec3::squared_norm(&v0);
    let l1 = vec3::squared_norm(&v1);
    let l2 = vec3::squared_norm(&v2);
    [
        (l1 + l2 - l0) * tmp, // cot0 = cos0/sin0 = {(l1*l1+l2*l2-l0*l0)/(2*l1*l2)} / {2*area/(l1*l2)}
        (l2 + l0 - l1) * tmp, // cot1 = cos1/sin1 = {(l2*l2+l0*l0-l1*l1)/(2*l2*l0)} / {2*area/(l2*l0)}
        (l0 + l1 - l2) * tmp, // cot2 = cos2/sin2 = {(l0*l0+l1*l1-l2*l2)/(2*l0*l1)} / {2*area/(l0*l1)}
    ]
}

pub fn emat_cotangent_laplacian<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [[[T; 1]; 3]; 3]
where
    T: num_traits::Float,
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
    T: num_traits::Float,
{
    use crate::vec3;
    let eps: T = T::epsilon();
    let edge1 = vec3::sub(p1, p0);
    let edge2 = vec3::sub(p2, p0);
    let pvec = vec3::cross(ray_dir, &edge2);
    let det = vec3::dot(&edge1, &pvec);
    if det > -eps && det < eps {
        return None;
    }
    let invdet = T::one() / det;
    let tvec = vec3::sub(ray_org, p0);
    let u = invdet * vec3::dot(&tvec, &pvec);
    if u < T::zero() || u > T::one() {
        return None;
    }
    let qvec = vec3::cross(&tvec, &edge1);
    let v = invdet * vec3::dot(ray_dir, &qvec);
    if v < T::zero() || u + v > T::one() {
        return None;
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    let t = invdet * vec3::dot(&edge2, &qvec);
    Some(t)
}

pub fn to_barycentric_coords<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3], q: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
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
    Real: num_traits::Float,
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

    pub fn normal(&self) -> [Real; 3] {
        normal(self.p0, self.p1, self.p2)
    }

    pub fn unit_normal(&self) -> [Real; 3] {
        let n = normal(self.p0, self.p1, self.p2);
        crate::vec3::normalized(&n)
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
