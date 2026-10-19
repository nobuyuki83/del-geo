//! Analytically integrated, unregularized 3D Kelvinlet on a flat triangle.
//!
//! Density is scalar, interpolated linearly from its three vertex values.
//! Force per unit area is density(y) * f; displacement is kernel(x) * f.
//! No normalization by total density or triangle area is performed.
//! Matrices are row-major: matrix[row][column].
//!
//! Supports points on the triangle, including its edges and vertices, by
//! analytic limits. This does NOT imply finite spatial derivatives of kernel.
//! Far-field cancellation and nearly degenerate triangles can reduce accuracy.
//!
//! Test: cargo test -p del-geo-core tri3_kelvinlet

use crate::vec3::{add, cross, dot, norm, scale, sub};

pub type Vec3<T> = [T; 3];
pub type Mat3ArrayOfRows<T> = [[T; 3]; 3];

#[derive(Debug, Clone, Copy)]
pub struct Evaluation<T> {
    /// P(x) = integral_T density(y) / |x-y| dA.
    pub potential: T,
    /// Hessian of S(x) = integral_T density(y) * |x-y| dA.
    pub hessian_s: Mat3ArrayOfRows<T>,
    /// Integrated elastic Green tensor; u(x) = kernel * f.
    pub kernel: Mat3ArrayOfRows<T>,
}

/// Evaluate the triangle Kelvinlet and Hessian of its scalar potential.
///
/// `mu` is the shear modulus (> 0); `nu` is Poisson's ratio (-1 < nu <= 0.5).
/// `nu == 0.5` evaluates the incompressible limit. Signed densities are allowed.
/// Input coordinates and density must be finite. Degenerate or extremely thin
/// triangles are rejected. Geometry is normalized internally by its longest edge.
pub fn evaluate<T>(
    vertices: [Vec3<T>; 3],
    density: [T; 3],
    x: Vec3<T>,
    mu: T,
    nu: T,
) -> Result<Evaluation<T>, &'static str>
where
    T: num_traits::Float + num_traits::FloatConst,
{
    if !vertices
        .iter()
        .flatten()
        .chain(x.iter())
        .chain(density.iter())
        .all(|a| a.is_finite())
    {
        return Err("coordinates and densities must be finite");
    }
    let zero = T::zero();
    let one = T::one();
    let two = one + one;
    let four = two + two;
    let half = one / two;
    let eight = two * two * two;
    let sixteen = eight + eight;
    let sixty_four = eight * eight;
    //
    if !mu.is_finite() || mu <= zero || !nu.is_finite() || nu <= -one || nu > half {
        return Err("require mu > 0 and -1 < nu <= 0.5");
    }
    let a = sub(&vertices[1], &vertices[0]);
    let b = sub(&vertices[2], &vertices[0]);
    let len_scale = norm(&a).max(norm(&b)).max(norm(&sub(&b, &a)));
    if !len_scale.is_finite() || len_scale == zero {
        return Err("invalid triangle scale");
    }
    let v = [
        [zero; 3],
        scale(&a, one / len_scale),
        scale(&b, one / len_scale),
    ];
    let xr = scale(&sub(&x, &vertices[0]), one / len_scale);
    if !xr.iter().all(|z| z.is_finite()) {
        return Err("normalized coordinates overflowed");
    }
    let normal = cross(&v[1], &v[2]);
    let area2 = norm(&normal);
    if area2 <= sixty_four * T::epsilon() {
        return Err("degenerate or numerically too thin triangle");
    }
    let n = scale(&normal, one / area2);
    let g = scale(
        &add(
            &scale(&cross(&v[2], &n), density[1] - density[0]),
            &scale(&cross(&n, &v[1]), density[2] - density[0]),
        ),
        one / area2,
    );

    let mut h = dot(&xr, &n);
    // Snap only roundoff-sized plane distances. This is not blob regularization.
    if h.abs() <= eight * T::epsilon() * (one + norm(&xr)) {
        h = zero;
    }
    let p = sub(&xr, &scale(&n, h));
    let c = density[0] + dot(&g, &p);

    let mut sum_dl = zero;
    let mut j_vec = [zero; 3];
    let mut h_d = [zero; 3];
    let mut h2_b = [zero; 3];
    let mut m_mat = [[zero; 3]; 3];
    for k in 0..3 {
        let va = v[k];
        let vb = v[(k + 1) % 3];
        let edge = sub(&vb, &va);
        let e = scale(&edge, one / norm(&edge));
        let m = cross(&e, &n);
        let ap = sub(&va, &p);
        let bp = sub(&vb, &p);
        let d = dot(&ap, &m);
        let t0 = dot(&ap, &e);
        let t1 = dot(&bp, &e);
        let q = h.hypot(d);
        let q2 = q * q;
        let r0 = t0.hypot(q);
        let r1 = t1.hypot(q);
        let (l, q2_l) = if q == zero {
            // L may diverge, but all its coefficients vanish here.
            // Store the limiting products d*L = q^2*L = h*L = 0.
            (zero, zero)
        } else {
            let l = asinh_scaled(t1, q) - asinh_scaled(t0, q);
            (l, q2 * l)
        };
        let a_edge = r1 - r0;
        let e_edge = half * (t1 * r1 - t0 * r0 + q2_l);
        let ce = c + d * dot(&g, &m);
        let ge = dot(&g, &e);
        let z = add(
            &scale(&e, ce * a_edge + ge * (e_edge - q2_l)),
            &scale(&m, d * (ce * l + ge * a_edge)),
        );
        sum_dl = sum_dl + d * l;
        j_vec = add(&j_vec, &scale(&m, e_edge));
        // Accumulate only products needed in the final tensor: no infinity*0.
        h_d = add(&h_d, &scale(&m, h * (ce * l + ge * a_edge)));
        h2_b = add(&h2_b, &scale(&m, h * h * l));
        for i in 0..3 {
            for j in 0..3 {
                m_mat[i][j] = m_mat[i][j] + z[i] * m[j];
            }
        }
    }

    let omega = if h == zero {
        // Only h*omega enters the formulas; its planar limit is zero.
        zero
    } else {
        let r = [sub(&xr, &v[0]), sub(&xr, &v[1]), sub(&xr, &v[2])];
        let len = [norm(&r[0]), norm(&r[1]), norm(&r[2])];
        let numerator = h * area2;
        let denominator = len[0] * len[1] * len[2]
            + dot(&r[0], &r[1]) * len[2]
            + dot(&r[1], &r[2]) * len[0]
            + dot(&r[2], &r[0]) * len[1];
        two * numerator.atan2(denominator)
    };
    let p0 = sum_dl - h * omega;
    let potential = c * p0 + dot(&g, &j_vec);
    let h_v = sub(&scale(&g, h * p0), &h_d);
    let w = c * h * omega - dot(&g, &h2_b);
    let mut hessian_s = [[zero; 3]; 3];
    let mut kernel = [[zero; 3]; 3];
    let coefficient = (one / mu) / (sixteen * T::PI() * (one - nu));
    for i in 0..3 {
        for j in 0..3 {
            let sym = half * (m_mat[i][j] + m_mat[j][i] - j_vec[i] * g[j] - g[i] * j_vec[j]);
            let hij = (potential - w) * n[i] * n[j] + sym + n[i] * h_v[j] + h_v[i] * n[j];
            hessian_s[i][j] = len_scale * hij;
            let diagonal = if i == j {
                four * (one - nu) * potential
            } else {
                zero
            };
            kernel[i][j] = len_scale * coefficient * (diagonal - hij);
        }
    }
    let out = Evaluation {
        potential: len_scale * potential,
        hessian_s,
        kernel,
    };
    if !out.potential.is_finite()
        || !out
            .hessian_s
            .iter()
            .flatten()
            .chain(out.kernel.iter().flatten())
            .all(|a| a.is_finite())
    {
        return Err("floating-point overflow; rescale inputs");
    }
    Ok(out)
}

/// Apply a row-major tensor to a force vector.
pub fn displacement<T>(kernel: Mat3ArrayOfRows<T>, force: Vec3<T>) -> Vec3<T>
where
    T: num_traits::Float,
{
    let kernel: &[T; 9] = kernel.as_flattened().try_into().unwrap();
    crate::mat3_row_major::mult_vec(kernel, &force)
}

fn asinh_scaled<T>(t: T, q: T) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    let ratio = t / q;
    if ratio.is_finite() {
        ratio.asinh()
    } else {
        t.signum() * (t.abs().ln() - q.ln() + T::LN_2())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type T = f64;
    const V: [Vec3<f64>; 3] = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.2, 0.9, 0.0]];
    fn close(a: f64, b: f64, tol: f64) {
        assert!(
            (a - b).abs() <= tol * (1.0 + a.abs().max(b.abs())),
            "{a} != {b}"
        );
    }
    // Independent quadrature, with Duffy coordinates about vertex 0.
    // The factor u cancels the point kernel's singularity at that vertex.
    fn quadrature(
        v: [Vec3<T>; 3],
        rho: Vec3<T>,
        x: Vec3<T>,
        steps: usize,
    ) -> (f64, Mat3ArrayOfRows<T>) {
        let a = sub(&v[1], &v[0]);
        let b = sub(&v[2], &v[0]);
        let jac = norm(&cross(&a, &b));
        let mut p = 0.0;
        let mut h = [[0.0; 3]; 3];
        for i in 0..steps {
            let u = (i as f64 + 0.5) / steps as f64;
            for j in 0..steps {
                let t = (j as f64 + 0.5) / steps as f64;
                let y = add(&v[0], &scale(&add(&scale(&a, 1.0 - t), &scale(&b, t)), u));
                let density = (1.0 - u) * rho[0] + u * ((1.0 - t) * rho[1] + t * rho[2]);
                let r = sub(&x, &y);
                let len = norm(&r);
                let weight = density * jac * u / (steps as f64).powi(2);
                p += weight / len;
                for row in 0..3 {
                    for col in 0..3 {
                        let id = if row == col { 1.0 } else { 0.0 };
                        h[row][col] += weight * (id / len - r[row] * r[col] / len.powi(3));
                    }
                }
            }
        }
        (p, h)
    }
    #[test]
    fn off_surface_matches_independent_quadrature() {
        let v = [[0.1, 0.2, 0.3], [1.3, 0.1, 0.6], [0.3, 1.1, 0.5]];
        for x in [[0.4, 0.25, 0.9], [-0.4, 0.7, -0.6]] {
            let out = evaluate::<T>(v, [0.8, 1.2, 0.5], x, 2.0, 0.3).unwrap();
            let (p, h) = quadrature(v, [0.8, 1.2, 0.5], x, 600);
            close(out.potential, p, 2e-6);
            for i in 0..3 {
                for j in 0..3 {
                    close(out.hessian_s[i][j], h[i][j], 2e-6);
                }
            }
        }
    }
    #[test]
    fn all_vertices_match_singular_quadrature() {
        let rho = [0.8, 1.2, 0.5];
        for k in 0..3 {
            let v = [V[k], V[(k + 1) % 3], V[(k + 2) % 3]];
            let r = [rho[k], rho[(k + 1) % 3], rho[(k + 2) % 3]];
            let out = evaluate::<T>(v, r, v[0], 1.0, 0.5).unwrap();
            let (p, h) = quadrature(v, r, v[0], 800);
            close(out.potential, p, 2e-6);
            for i in 0..3 {
                for j in 0..3 {
                    close(out.hessian_s[i][j], h[i][j], 2e-6);
                }
            }
        }
    }
    #[test]
    fn plane_limits_and_winding() {
        for xy in [[0.3, 0.2, 0.0], [0.6, 0.0, 0.0], V[0], V[1], V[2]] {
            let at = evaluate::<T>(V, [1.0; 3], xy, 1.0, 0.3).unwrap();
            let reversed = evaluate::<T>([V[0], V[2], V[1]], [1.0; 3], xy, 1.0, 0.3).unwrap();
            for sign in [-1.0, 1.0] {
                let near = evaluate::<T>(V, [1.0; 3], add(&xy, &[0.0, 0.0, sign * 1e-9]), 1.0, 0.3)
                    .unwrap();
                for i in 0..3 {
                    for j in 0..3 {
                        close(at.kernel[i][j], reversed.kernel[i][j], 1e-12);
                        close(at.kernel[i][j], near.kernel[i][j], 2e-7);
                    }
                }
            }
        }
    }
    #[test]
    fn constant_density_vertex_has_known_potential() {
        let v = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let out = evaluate::<T>(v, [1.0; 3], v[0], 1.0, 0.3).unwrap();
        close(out.potential, 2.0_f64.sqrt() * 1.0_f64.asinh(), 1e-13);
        close(
            (0..3).map(|i| out.hessian_s[i][i]).sum(),
            2.0 * out.potential,
            1e-13,
        );
    }
    #[test]
    fn scale_density_and_invalid_inputs() {
        let x = [0.4, 0.2, 0.7];
        let base = evaluate::<T>(V, [0.8, 1.2, 0.5], x, 1.0, 0.3).unwrap();
        let scaled = evaluate::<T>(
            V.map(|v| scale(&v, 7.0)),
            [1.6, 2.4, 1.0],
            scale(&x, 7.0),
            1.0,
            0.3,
        )
        .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(scaled.kernel[i][j], 14.0 * base.kernel[i][j], 1e-12);
            }
        }
        assert!(evaluate::<T>([[0.0; 3]; 3], [1.0; 3], x, 1.0, 0.3).is_err());
        assert!(evaluate::<T>(V, [1.0; 3], x, 0.0, 0.3).is_err());
        assert!(evaluate::<T>(V, [1.0; 3], x, 1.0, 0.6).is_err());
    }
}
