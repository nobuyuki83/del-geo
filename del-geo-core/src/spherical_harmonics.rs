#![allow(clippy::needless_range_loop)]
#![allow(clippy::comparison_chain)]

use num_complex::{Complex, ComplexFloat};
use std::f64::consts::PI;

/// Calculate the normalization of the vector.
#[inline]
pub fn normalize(x: &mut f64, y: &mut f64, z: &mut f64) -> f64 {
    let r = f64::sqrt(*x * *x + *y * *y + *z * *z);
    let invr = 1.0 / r;
    *x *= invr;
    *y *= invr;
    *z *= invr;
    r
}

/// Calculate the coefficients of the spherical harmonics for l <= 9 and store them in an array buffer.
/// Try to access the coefficient Y_l^m by the index: base + l + m, where base = l^2.
#[inline]
pub fn sph_coeff_buffer(n: i8, x: f64, y: f64, z: f64) -> [f64; 100] {
    let inv_pi = 1.0 / PI;
    let ep = Complex::new(x, y);
    let mut res = [0.0; 100];

    // n = 0, single spheric, s orbital
    res[0] = 0.5 * f64::sqrt(inv_pi);
    if n == 0 {
        return res;
    }

    // n = 1, spindle, p orbital
    let r1 = ep.re();
    let i1 = ep.im();
    let v1 = -0.5 * f64::sqrt(1.5 * inv_pi);
    res[1] = v1 * i1;
    res[2] = 0.5 * f64::sqrt(3.0 * inv_pi) * z;
    res[3] = v1 * r1;
    if n == 1 {
        return res;
    }

    // n = 2, double spindle, d orbital
    let ep2 = ep * ep;
    let r2 = ep2.re();
    let i2 = ep2.im();
    let v1 = -0.5 * f64::sqrt(7.5 * inv_pi) * z;
    let v2 = 0.25 * f64::sqrt(7.5 * inv_pi);
    res[4] = v2 * i2;
    res[5] = v1 * i1;
    res[6] = 0.25 * f64::sqrt(5.0 * inv_pi) * (2.0 * z * z - x * x - y * y);
    res[7] = v1 * r1;
    res[8] = v2 * r2;
    if n == 2 {
        return res;
    }

    // n = 3, triple spindle, f orbital
    let ep3 = ep2 * ep;
    let r3 = ep3.re();
    let i3 = ep3.im();
    let v1 = -0.125 * f64::sqrt(21.0 * inv_pi) * (4.0 * z * z - x * x - y * y);
    let v2 = 0.250 * f64::sqrt(52.5 * inv_pi) * z;
    let v3 = -0.125 * f64::sqrt(35.0 * inv_pi);
    res[9] = v3 * i3;
    res[10] = v2 * i2;
    res[11] = v1 * i1;
    res[12] = 0.250 * f64::sqrt(7.0 * inv_pi) * z * (-3.0 * x * x - 3.0 * y * y + 2.0 * z * z);
    res[13] = v1 * r1;
    res[14] = v2 * r2;
    res[15] = v3 * r3;
    if n == 3 {
        return res;
    }

    // n = 4, quadruple spindle, g orbital
    let ep4 = ep3 * ep;
    let r4 = ep4.re();
    let i4 = ep4.im();
    let z2 = z * z;
    let v1 = -3.0 / 8.00 * f64::sqrt(5.0 * inv_pi) * z * (7.0 * z2 - 3.0);
    let v2 = 3.0 / 8.00 * f64::sqrt(5.0 * 0.5 * inv_pi) * (7.0 * z2 - 1.0);
    let v3 = -3.0 / 8.00 * f64::sqrt(35.0 * inv_pi) * z;
    let v4 = 3.0 / 16.0 * f64::sqrt(35.0 * 0.5 * inv_pi);
    res[16] = v4 * i4;
    res[17] = v3 * i3;
    res[18] = v2 * i2;
    res[19] = v1 * i1;
    res[20] = 3.0 / 16.0 * f64::sqrt(inv_pi) * (35.0 * z2 * z2 - 30.0 * z2 + 3.0);
    res[21] = v1 * r1;
    res[22] = v2 * r2;
    res[23] = v3 * r3;
    res[24] = v4 * r4;
    if n == 4 {
        return res;
    }

    // n = 5, quintuple spindle, h orbital
    let r5 = (ep4 * ep).re();
    let i5 = (ep4 * ep).im();
    let z4 = z2 * z2;
    let v1 = -1.0 / 16.0 * f64::sqrt(82.5 * inv_pi) * (21.0 * z4 - 14.0 * z2 + 1.0);
    let v2 = 1.0 / 8.00 * f64::sqrt(577.5 * inv_pi) * z * (3.0 * z2 - 1.0);
    let v3 = -1.0 / 32.0 * f64::sqrt(385.0 * inv_pi) * (9.0 * z2 - 1.0);
    let v4 = 3.0 / 16.0 * f64::sqrt(192.5 * inv_pi) * z;
    let v5 = -3.0 / 32.0 * f64::sqrt(77.0 * inv_pi);
    res[25] = v5 * i5;
    res[26] = v4 * i4;
    res[27] = v3 * i3;
    res[28] = v2 * i2;
    res[29] = v1 * i1;
    res[30] = 1.0 / 16.0 * f64::sqrt(11.0 * inv_pi) * z * (63.0 * z4 - 70.0 * z2 + 15.0);
    res[31] = v1 * r1;
    res[32] = v2 * r2;
    res[33] = v3 * r3;
    res[34] = v4 * r4;
    res[35] = v5 * r5;
    if n == 5 {
        return res;
    }

    // n = 6, sextuple spindle, i orbital
    let r6 = (ep4 * ep2).re();
    let i6 = (ep4 * ep2).im();
    let v1 = -1.0 / 16.0 * f64::sqrt(273.0 * 0.5 * inv_pi) * z * (33.0 * z4 - 30.0 * z2 + 5.0);
    let v2 = 1.0 / 64.0 * f64::sqrt(1365.0 * inv_pi) * (33.0 * z4 - 18.0 * z2 + 1.0);
    let v3 = -1.0 / 32.0 * f64::sqrt(1365.0 * inv_pi) * z * (11.0 * z2 - 3.0);
    let v4 = 3.0 / 32.0 * f64::sqrt(91.0 * 0.5 * inv_pi) * (11.0 * z2 - 1.0);
    let v5 = -3.0 / 32.0 * f64::sqrt(1001.0 * inv_pi) * z;
    let v6 = 1.0 / 64.0 * f64::sqrt(3003.0 * inv_pi);
    res[36] = v6 * i6;
    res[37] = v5 * i5;
    res[38] = v4 * i4;
    res[39] = v3 * i3;
    res[40] = v2 * i2;
    res[41] = v1 * i1;
    res[42] =
        1.0 / 32.0 * f64::sqrt(13.0 * inv_pi) * (231.0 * z4 * z2 - 315.0 * z4 + 105.0 * z2 - 5.0);
    res[43] = v1 * r1;
    res[44] = v2 * r2;
    res[45] = v3 * r3;
    res[46] = v4 * r4;
    res[47] = v5 * r5;
    res[48] = v6 * r6;
    if n == 6 {
        return res;
    }
    // n = 7, septuple spindle, j orbital
    let r7 = (ep4 * ep2 * ep).re();
    let i7 = (ep4 * ep2 * ep).im();
    let v1 = -1.0 / 64.0
        * f64::sqrt(105.0 * 0.5 * inv_pi)
        * (429.0 * z4 * z2 - 495.0 * z4 + 135.0 * z2 - 5.0);
    let v2 = 3.0 / 64.0 * f64::sqrt(35.0 * inv_pi) * (143.0 * z4 * z - 110.0 * z2 * z + 15.0 * z);
    let v3 = -3.0 / 64.0 * f64::sqrt(35.0 * 0.5 * inv_pi) * (143.0 * z4 - 66.0 * z2 + 3.0);
    let v4 = 3.0 / 32.0 * f64::sqrt(385.0 * 0.5 * inv_pi) * (13.0 * z2 * z - 3.0 * z);
    let v5 = -3.0 / 64.0 * f64::sqrt(385.0 * 0.5 * inv_pi) * (13.0 * z2 - 1.0);
    let v6 = 3.0 / 64.0 * f64::sqrt(5005.0 * inv_pi) * z;
    let v7 = -3.0 / 128.0 * f64::sqrt(1430.0 * inv_pi);
    res[49] = v7 * i7;
    res[50] = v6 * i6;
    res[51] = v5 * i5;
    res[52] = v4 * i4;
    res[53] = v3 * i3;
    res[54] = v2 * i2;
    res[55] = v1 * i1;
    res[56] = 1.0 / 32.0
        * f64::sqrt(15.0 * inv_pi)
        * (429.0 * z4 * z2 * z - 693.0 * z4 * z + 315.0 * z2 * z - 35.0 * z);
    res[57] = v1 * r1;
    res[58] = v2 * r2;
    res[59] = v3 * r3;
    res[60] = v4 * r4;
    res[61] = v5 * r5;
    res[62] = v6 * r6;
    res[63] = v7 * r7;
    if n == 7 {
        return res;
    }
    // n = 8, octuple spindle, k orbital
    let ep8 = ep4 * ep4;
    let r8 = ep8.re();
    let i8 = ep8.im();
    let z8 = z4 * z4;
    let v1 = -3.0 / 64.00
        * f64::sqrt(17.0 * 0.5 * inv_pi)
        * (715.0 * z4 * z2 * z - 1001.0 * z4 * z + 385.0 * z2 * z - 35.0 * z);
    let v2 =
        3.0 / 128.0 * f64::sqrt(595.0 * inv_pi) * (143.0 * z4 * z2 - 143.0 * z4 + 33.0 * z2 - 1.0);
    let v3 = -1.0 / 64.00
        * f64::sqrt(19635.0 * 0.5 * inv_pi)
        * (39.0 * z4 * z - 26.0 * z2 * z + 3.0 * z);
    let v4 = 3.0 / 128.0 * f64::sqrt(1309.0 * 0.5 * inv_pi) * (65.0 * z4 - 26.0 * z2 + 1.0);
    let v5 = -3.0 / 64.00 * f64::sqrt(17017.0 * 0.5 * inv_pi) * (5.0 * z2 * z - z);
    let v6 = 1.0 / 128.0 * f64::sqrt(7293.0 * inv_pi) * (15.0 * z2 - 1.0);
    let v7 = -3.0 / 64.00 * f64::sqrt(12155.0 * 0.5 * inv_pi) * z;
    let v8 = 3.0 / 256.0 * f64::sqrt(12155.0 * 0.5 * inv_pi);
    res[64] = v8 * i8;
    res[65] = v7 * i7;
    res[66] = v6 * i6;
    res[67] = v5 * i5;
    res[68] = v4 * i4;
    res[69] = v3 * i3;
    res[70] = v2 * i2;
    res[71] = v1 * i1;
    res[72] = 1.0 / 256.0
        * f64::sqrt(17.0 * inv_pi)
        * (6435.0 * z8 - 12012.0 * z4 * z2 + 6930.0 * z4 - 1260.0 * z2 + 35.0);
    res[73] = v1 * r1;
    res[74] = v2 * r2;
    res[75] = v3 * r3;
    res[76] = v4 * r4;
    res[77] = v5 * r5;
    res[78] = v6 * r6;
    res[79] = v7 * r7;
    res[80] = v8 * r8;
    if n == 8 {
        return res;
    }
    // n = 9, nonuple spindle, l orbital
    let r9 = (ep8 * ep).re();
    let i9 = (ep8 * ep).im();
    let v1 = -3.0 / 256.0
        * f64::sqrt(95.0 * 0.5 * inv_pi)
        * (2431.0 * z8 - 4004.0 * z4 * z2 + 2002.0 * z4 - 308.0 * z2 + 7.0);
    let v2 = 3.0 / 128.0
        * f64::sqrt(1045.0 * inv_pi)
        * z
        * (221.0 * z4 * z2 - 273.0 * z4 + 91.0 * z2 - 7.0);
    let v3 = -1.0 / 256.0
        * f64::sqrt(21945.0 * inv_pi)
        * (221.0 * z4 * z2 - 195.0 * z4 + 39.0 * z2 - 1.0);
    let v4 = 3.0 / 256.0 * f64::sqrt(95095.0 * 2.0 * inv_pi) * z * (17.0 * z4 - 10.0 * z2 + 1.0);
    let v5 = -3.0 / 256.0 * f64::sqrt(2717.0 * inv_pi) * (85.0 * z4 - 30.0 * z2 + 1.0);
    let v6 = 1.0 / 128.0 * f64::sqrt(40755.0 * inv_pi) * z * (17.0 * z2 - 3.0);
    let v7 = -3.0 / 512.0 * f64::sqrt(13585.0 * inv_pi) * (17.0 * z2 - 1.0);
    let v8 = 3.0 / 512.0 * f64::sqrt(230945.0 * 2.0 * inv_pi) * z;
    let v9 = -1.0 / 512.0 * f64::sqrt(230945.0 * inv_pi);
    res[81] = v9 * i9;
    res[82] = v8 * i8;
    res[83] = v7 * i7;
    res[84] = v6 * i6;
    res[85] = v5 * i5;
    res[86] = v4 * i4;
    res[87] = v3 * i3;
    res[88] = v2 * i2;
    res[89] = v1 * i1;
    res[90] = 1.0 / 256.0
        * f64::sqrt(19.0 * inv_pi)
        * z
        * (12155.0 * z8 - 25740.0 * z4 * z2 + 18018.0 * z4 - 4620.0 * z2 + 315.0);
    res[91] = v1 * r1;
    res[92] = v2 * r2;
    res[93] = v3 * r3;
    res[94] = v4 * r4;
    res[95] = v5 * r5;
    res[96] = v6 * r6;
    res[97] = v7 * r7;
    res[98] = v8 * r8;
    res[99] = v9 * r9;
    res
}

/// Calculate the coefficients of Legendre Polynomials in all orders <= n.
/// Try to access the coefficient of x^m in P_l(x) by result[l][m].
fn legendre_coeff_vec(n: u64) -> Vec<Vec<f64>> {
    let mut res = Vec::new();

    let p0 = vec![1.0];
    res.push(p0);

    let p1 = vec![0.0, 1.0];
    res.push(p1);

    for i in 2..n + 1 {
        let mut p: Vec<f64> = Vec::new();
        for j in 0..i + 1 {
            if (i - j) % 2 != 0 {
                p.push(0.0);
            } else if j == 0 {
                let inv_i = 1.0 / i as f64;
                p.push(-(1.0 - inv_i) * res[i as usize - 2][j as usize]);
            } else {
                let inv_i = 1.0 / i as f64;
                let tmp = if i >= 1 && j >= 1 && j <= i {
                    (2.0 - inv_i) * res[i as usize - 1][j as usize - 1]
                } else {
                    0.0
                } - if i >= 2 && j <= i - 2 {
                    (1.0 - inv_i) * res[i as usize - 2][j as usize]
                } else {
                    0.0
                };
                p.push(tmp);
            }
        }
        res.push(p);
    }
    res
}

/// Calculate the factorial of a number.
/// Here needs a lot of optimization to avoid overflow. Now the maximum parameter is 11, which is enough for the current use.
fn factorial(n: u128) -> u128 {
    if n == 0 {
        return 1;
    }
    let mut res = 1;
    for i in 1..n + 1 {
        res *= i;
    }
    res
}

/// Calculate the coefficient of the term in the Legendre Polynomial.
/// The term is P_l^m(x) = get_legendre_poly_term_coeff(l, m) * x^m.
pub fn get_legendre_poly_term_coeff(func_order: u32, term_order: u32) -> f64 {
    if (func_order - term_order) % 2 != 0 {
        0.0
    } else {
        let k = (func_order - term_order) / 2;
        let mol = if k % 2 != 0 { -1.0 } else { 1.0 }
            * factorial((func_order + term_order) as u128) as f64;
        let den = (2_u32.pow(func_order) as u128
            * factorial(k as u128)
            * factorial(term_order as u128)
            * factorial((k + term_order).into())) as f64;
        mol / den
    }
}

/// Calculate the associated Legendre Polynomial P_l^m(x) for l >= |m|.
fn calculate_assoc_legendre_poly(l: u64, m: i64, x: f64) -> f64 {
    let m_abs = m.abs();
    assert!(l >= m_abs as u64);
    assert!((-1.0..=1.0).contains(&x));
    let sign = m.signum();
    let legendre = legendre_coeff_vec(l);
    let legendre_coeff_vec = &legendre[l as usize];
    let mut sum_coeff: Vec<f64> = Vec::new();
    for i in m_abs as usize..l as usize + 1 {
        let mut tmp = legendre_coeff_vec[i];
        tmp *= factorial(i as u128) as f64 / factorial((i as u32 - m_abs as u32).into()) as f64;
        sum_coeff.push(tmp);
    }
    let mut sum = 0.0;
    for i in 0..sum_coeff.len() {
        sum += sum_coeff[i] * f64::powi(x, i as i32);
    }
    let mut res =
        if m_abs % 2 == 0 { 1.0 } else { -1.0 } * (1.0 - x * x).powf(m_abs as f64 / 2.0) * sum;
    res *= if sign == 1 {
        1.0
    } else {
        (if m.abs() % 2 == 0 { 1.0 } else { -1.0 })
            * factorial((l - m.unsigned_abs()) as u128) as f64
            / factorial((l + m.unsigned_abs()) as u128) as f64
    } as f64;
    res
}

/// Calculate the coefficient of the spherical harmonics for "any" l and m.
/// However, u128 is not enough for the factorial calculation, so the maximum l is 11.
/// Now that it is related to multiply and divide in big factorials, it has large opmitization potential.
/// But optimization based on combination is hard in both mathematically and programmatically.
pub fn get_spherical_harmonics_coeff(l: i64, m: i64, x: f64, y: f64, z: f64) -> f64 {
    let m_abs = m.abs();
    assert!(l >= 0 && l >= m_abs);
    let r = f64::sqrt(x * x + y * y);
    let ep = Complex::new(x / r, y / r);
    if m == 0 {
        f64::sqrt((2.0 * l as f64 + 1.0) / 4.0 / PI) * calculate_assoc_legendre_poly(l as u64, 0, z)
    } else if m < 0 {
        f64::sqrt(
            (2.0 * l as f64 + 1.0) / (4.0 * PI)
                * factorial((l as u32 - m_abs as u32).into()) as f64
                / factorial((l as u32 + m_abs as u32).into()) as f64,
        ) * calculate_assoc_legendre_poly(l as u64, m_abs, z)
            * ep.powf(m_abs as f64).im()
    } else {
        f64::sqrt(
            (2.0 * l as f64 + 1.0) * factorial((l as u32 - m_abs as u32).into()) as f64
                / factorial((l as u32 + m_abs as u32).into()) as f64
                / (4.0 * PI),
        ) * calculate_assoc_legendre_poly(l as u64, m_abs, z)
            * ep.powf(m_abs as f64).re()
    }
}

#[test]
fn test_get_spherical_harmonics_coeff() {
    // let tmp = 1.0 / f64::sqrt(3.0);
    let sph_coeff = sph_coeff_buffer(
        9,
        1.0 / f64::sqrt(3.0),
        1.0 / f64::sqrt(3.0),
        1.0 / f64::sqrt(3.0),
    );
    for l in 7..10 {
        for m in -l..l + 1 {
            let current_coeff = get_spherical_harmonics_coeff(
                l,
                m,
                1.0 / f64::sqrt(3.0),
                1.0 / f64::sqrt(3.0),
                1.0 / f64::sqrt(3.0),
            );
            let base = l.pow(2);
            assert!(base + m >= 0 && base + l + m < 100);
            assert!((current_coeff - sph_coeff[(base + l + m) as usize]).abs() <= 1.0e-7);
        }
    }
}
