//! methods for hexahedron.
//! the coordinates of the points are stored in the array of array as`[[Real;3];8]`.
//! where `[xyz, Xyz, XYz, xYz, xyZ, XyZ, XYZ, xYZ]`.

pub fn shapefunc<Real>(
    node2xyz: &[[Real; 3]; 8],
    r0: Real,
    r1: Real,
    r2: Real,
) -> ([Real; 8], [[Real; 3]; 8], Real)
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let one8 = one / (two * two * two);
    let an = [
        one8 * (one - r0) * (one - r1) * (one - r2), // xyz
        one8 * (one + r0) * (one - r1) * (one - r2), // Xyz
        one8 * (one + r0) * (one + r1) * (one - r2), // XYz
        one8 * (one - r0) * (one + r1) * (one - r2), // xYz
        one8 * (one - r0) * (one - r1) * (one + r2),
        one8 * (one + r0) * (one - r1) * (one + r2),
        one8 * (one + r0) * (one + r1) * (one + r2),
        one8 * (one - r0) * (one + r1) * (one + r2),
    ];

    let dndr = [
        [
            -one8 * (one - r1) * (one - r2),
            -one8 * (one - r0) * (one - r2),
            -one8 * (one - r0) * (one - r1),
        ],
        [
            one8 * (one - r1) * (one - r2),
            -one8 * (one + r0) * (one - r2),
            -one8 * (one + r0) * (one - r1),
        ],
        [
            one8 * (one + r1) * (one - r2),
            one8 * (one + r0) * (one - r2),
            -one8 * (one + r0) * (one + r1),
        ],
        [
            -one8 * (one + r1) * (one - r2),
            one8 * (one - r0) * (one - r2),
            -one8 * (one - r0) * (one + r1),
        ],
        [
            -one8 * (one - r1) * (one + r2),
            -one8 * (one - r0) * (one + r2),
            one8 * (one - r0) * (one - r1),
        ],
        [
            one8 * (one - r1) * (one + r2),
            -one8 * (one + r0) * (one + r2),
            one8 * (one + r0) * (one - r1),
        ],
        [
            one8 * (one + r1) * (one + r2),
            one8 * (one + r0) * (one + r2),
            one8 * (one + r0) * (one + r1),
        ],
        [
            -one8 * (one + r1) * (one + r2),
            one8 * (one - r0) * (one + r2),
            one8 * (one - r0) * (one + r1),
        ],
    ];

    let (dndx, detjac) = crate::hex::grad_shapefunc_from_dndr(node2xyz, &dndr);
    (an, dndx, detjac)
}

pub fn grad_shapefunc_from_dndr<Real>(
    node2xyz: &[[Real; 3]; 8],
    dndr: &[[Real; 3]; 8],
) -> ([[Real; 3]; 8], Real)
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let mut dxdr = [[zero; 3]; 3];
    for inode in 0..8 {
        for idim in 0..3 {
            for jdim in 0..3 {
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[inode][idim] * dndr[inode][jdim];
            }
        }
    }

    let detjac = dxdr[0][0] * dxdr[1][1] * dxdr[2][2]
        + dxdr[1][0] * dxdr[2][1] * dxdr[0][2]
        + dxdr[2][0] * dxdr[0][1] * dxdr[1][2]
        - dxdr[0][0] * dxdr[2][1] * dxdr[1][2]
        - dxdr[1][0] * dxdr[0][1] * dxdr[2][2]
        - dxdr[2][0] * dxdr[1][1] * dxdr[0][2];

    let inv_jac = Real::one() / detjac;

    let drdx = [
        [
            inv_jac * (dxdr[1][1] * dxdr[2][2] - dxdr[1][2] * dxdr[2][1]),
            inv_jac * (dxdr[0][2] * dxdr[2][1] - dxdr[0][1] * dxdr[2][2]),
            inv_jac * (dxdr[0][1] * dxdr[1][2] - dxdr[0][2] * dxdr[1][1]),
        ],
        [
            inv_jac * (dxdr[1][2] * dxdr[2][0] - dxdr[1][0] * dxdr[2][2]),
            inv_jac * (dxdr[0][0] * dxdr[2][2] - dxdr[0][2] * dxdr[2][0]),
            inv_jac * (dxdr[0][2] * dxdr[1][0] - dxdr[0][0] * dxdr[1][2]),
        ],
        [
            inv_jac * (dxdr[1][0] * dxdr[2][1] - dxdr[1][1] * dxdr[2][0]),
            inv_jac * (dxdr[0][1] * dxdr[2][0] - dxdr[0][0] * dxdr[2][1]),
            inv_jac * (dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0]),
        ],
    ];

    let mut dndx = [[zero; 3]; 8];
    for inode in 0..8 {
        dndx[inode][0] =
            dndr[inode][0] * drdx[0][0] + dndr[inode][1] * drdx[1][0] + dndr[inode][2] * drdx[2][0];
        dndx[inode][1] =
            dndr[inode][0] * drdx[0][1] + dndr[inode][1] * drdx[1][1] + dndr[inode][2] * drdx[2][1];
        dndx[inode][2] =
            dndr[inode][0] * drdx[0][2] + dndr[inode][1] * drdx[1][2] + dndr[inode][2] * drdx[2][2];
    }

    (dndx, detjac)
}

pub fn grad_shapefunc<Real>(
    node2xyz: &[[Real; 3]; 8],
    quadrature: &[[Real; 2]],
    ir1: usize,
    ir2: usize,
    ir3: usize,
) -> ([[Real; 3]; 8], Real)
where
    Real: num_traits::Float,
{
    let r1 = quadrature[ir1][0];
    let r2 = quadrature[ir2][0];
    let r3 = quadrature[ir3][0];
    let (_an, dndx, detjac) = crate::hex::shapefunc(node2xyz, r1, r2, r3);
    let detwei = detjac * quadrature[ir1][1] * quadrature[ir2][1] * quadrature[ir3][1];
    (dndx, detwei)
}
