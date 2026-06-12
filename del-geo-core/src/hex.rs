//! methods for hexahedron.
//! the coordinates of the points are stored in the array of array as`[[Real;3];8]`.
//! where `[xyz, Xyz, XYz, xYz, xyZ, XyZ, XYZ, xYZ]`.

pub const HEX_SIGN: [[f64; 3]; 8] = [
    [-1., -1., -1.],
    [1., -1., -1.],
    [1., 1., -1.],
    [-1., 1., -1.],
    [-1., -1., 1.],
    [1., -1., 1.],
    [1., 1., 1.],
    [-1., 1., 1.],
];

pub const EDGE2NODE: [[usize; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

pub const FACE2IDX: [usize; 7] = [0, 4, 8, 12, 16, 20, 24];
pub const IDX2NODE: [usize; 24] = [
    0, 3, 2, 1, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7, 4, 5, 6, 7,
];

pub fn shapefunc<Real>(pco: &[Real; 3]) -> [Real; 8]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let one8 = one / (two * two * two);
    let (r0, r1, r2) = (pco[0], pco[1], pco[2]);
    [
        one8 * (one - r0) * (one - r1) * (one - r2), // xyz
        one8 * (one + r0) * (one - r1) * (one - r2), // Xyz
        one8 * (one + r0) * (one + r1) * (one - r2), // XYz
        one8 * (one - r0) * (one + r1) * (one - r2), // xYz
        one8 * (one - r0) * (one - r1) * (one + r2),
        one8 * (one + r0) * (one - r1) * (one + r2),
        one8 * (one + r0) * (one + r1) * (one + r2),
        one8 * (one - r0) * (one + r1) * (one + r2),
    ]
}

pub fn dndr<Real>(pco: &[Real; 3]) -> [[Real; 3]; 8]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let two = one + one;
    let one8 = one / (two * two * two);
    let (r0, r1, r2) = (pco[0], pco[1], pco[2]);
    [
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
    ]
}

pub fn dndx<Real>(
    node2xyz: &[[Real; 3]; 8],
    r0: Real,
    r1: Real,
    r2: Real,
) -> ([Real; 8], [[Real; 3]; 8], Real)
where
    Real: num_traits::Float,
{
    let pco = [r0, r1, r2];
    let dndr = dndr(&pco);
    let an = shapefunc(&pco);
    let (dndx, detjac) = crate::hex::grad_shapefunc_from_dndr(node2xyz, &dndr);
    (an, dndx, detjac)
}

pub fn dxdr<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    p6: &[Real; 3],
    p7: &[Real; 3],
    pco: &[Real; 3],
) -> [[Real; 3]; 3]
where
    Real: num_traits::Float,
{
    let dndr = dndr(pco);
    let zero = Real::zero();
    let mut dxdr = [[zero; 3]; 3];
    for idim in 0..3 {
        for jdim in 0..3 {
            dxdr[idim][jdim] = dxdr[idim][jdim] + p0[idim] * dndr[0][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p1[idim] * dndr[1][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p2[idim] * dndr[2][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p3[idim] * dndr[3][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p4[idim] * dndr[4][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p5[idim] * dndr[5][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p6[idim] * dndr[6][jdim];
            dxdr[idim][jdim] = dxdr[idim][jdim] + p7[idim] * dndr[7][jdim];
        }
    }
    dxdr
}

pub fn grad_shapefunc_from_dndr<Real>(
    node2xyz: &[[Real; 3]; 8],
    dndr: &[[Real; 3]; 8],
) -> ([[Real; 3]; 8], Real)
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let dxdr = {
        let mut dxdr = [[zero; 3]; 3];
        for idim in 0..3 {
            for jdim in 0..3 {
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[0][idim] * dndr[0][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[1][idim] * dndr[1][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[2][idim] * dndr[2][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[3][idim] * dndr[3][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[4][idim] * dndr[4][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[5][idim] * dndr[5][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[6][idim] * dndr[6][jdim];
                dxdr[idim][jdim] = dxdr[idim][jdim] + node2xyz[7][idim] * dndr[7][jdim];
            }
        }
        dxdr
    };

    let zero = Real::zero();
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
    let (_an, dndx, detjac) = crate::hex::dndx(node2xyz, r1, r2, r3);
    let detwei = detjac * quadrature[ir1][1] * quadrature[ir2][1] * quadrature[ir3][1];
    (dndx, detwei)
}

const ISO_EDGE: [i32; 256] = [
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03,
    0xe09, 0xf00, 0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f,
    0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6,
    0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569,
    0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69,
    0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6,
    0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c, 0xe5c,
    0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf,
    0x1c5, 0xcc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3,
    0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55, 0x35f, 0x256, 0x55a,
    0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5,
    0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65,
    0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa,
    0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0, 0xd30,
    0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33,
    0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f,
    0x596, 0x29a, 0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0,
];

const ISO_TRI: [[i8; 16]; 256] = [
    [
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    ],
    [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
    [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
    [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
    [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
    [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
    [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
    [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
    [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
    [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
    [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
    [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
    [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
    [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
    [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
    [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
    [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
    [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
    [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
    [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
    [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
    [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
    [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
    [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
    [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
    [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
    [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
    [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
    [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
    [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
    [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
    [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
    [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
    [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
    [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
    [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
    [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
    [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
    [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
    [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
    [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
    [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
    [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
    [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
    [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
    [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
    [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
    [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
    [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
    [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
    [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
    [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
    [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
    [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
    [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
    [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
    [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
    [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
    [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
    [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
    [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
    [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
    [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
    [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
    [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
    [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
    [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
    [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
    [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
    [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
    [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
    [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
    [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
    [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
    [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
    [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
    [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
    [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
    [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
    [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
    [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
    [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
    [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
    [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
    [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
    [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
    [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
    [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
    [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
    [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
    [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
    [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
    [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
    [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
    [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
    [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
    [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
    [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
    [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
    [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
    [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
    [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
    [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
    [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
    [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
    [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
    [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
    [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
    [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
    [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
    [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
    [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
    [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
    [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
    [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
    [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
    [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
    [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
    [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
    [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
    [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
    [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
    [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
    [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
    [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
    [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
    [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
    [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
    [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
    [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
    [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
    [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
    [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
    [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
    [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
    [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
    [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
    [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
    [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
    [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
    [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
    [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
    [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
    [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
    [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
    [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
    [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
    [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
    [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
    [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
    [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
    [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
    [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
    [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
    [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
    [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
    [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
    [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
    [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
    [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
    [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
    [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
    [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
    [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
    [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
    [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    ],
];

/// Interpolates the intersection point of an isosurface with an edge of a
/// hexahedral cell,
///
/// * `cent`  – center of the cell (x, y, z).
/// * `hw`    – half-width of the cell.
/// * `ind0` / `ind1` – indices of the two cube vertices that form the edge being inspected (0‥7, following `EDGE2NODE`).
/// * `dist0` / `dist1` – signed distances at the two end-points.
///
/// Returns the interpolated 3-D point as `[f64; 3]`.
pub fn iso_position_on_edge(
    pos_iso: &mut [f64; 3],
    cent: [f64; 3],
    hw: f64,
    i0_node: usize,
    i1_node: usize,
    val0: f64,
    val1: f64,
) {
    /// Vertex offsets for a unit cube, laid out in the same order as Marching-Cubes
    /// (x, y, z ∈ {0, 1}).  Adapt as needed if your ordering differs.
    use crate::vec3::Vec3;
    let p0 = HEX_SIGN[i0_node].scale(hw).add(&cent);
    let p1 = HEX_SIGN[i1_node].scale(hw).add(&cent);
    let denom = val1 - val0;
    let r0 = val1 / denom;
    let r1 = -val0 / denom;
    pos_iso[0] = p0[0] * r0 + p1[0] * r1;
    pos_iso[1] = p0[1] * r0 + p1[1] * r1;
    pos_iso[2] = p0[2] * r0 + p1[2] * r1;
}

pub fn volume<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    p6: &[Real; 3],
    p7: &[Real; 3],
    i_gauss_degree: usize,
) -> Real
where
    Real: num_traits::Float + 'static + std::fmt::Debug + std::iter::Sum,
    crate::quadrature_line::Quad<Real>: crate::quadrature_line::QuadratureLine<Real>,
{
    use crate::quadrature_line::QuadratureLine;
    let quadrature: &[[Real; 2]] = crate::quadrature_line::Quad::<Real>::hoge(i_gauss_degree);
    let nq = quadrature.len();
    let vol = itertools::iproduct!(0..nq, 0..nq, 0..nq)
        .map(|(ir0, ir1, ir2)| {
            let r0 = quadrature[ir0][0];
            let r1 = quadrature[ir1][0];
            let r2 = quadrature[ir2][0];
            let w = quadrature[ir0][1] * quadrature[ir1][1] * quadrature[ir2][1];
            let dxdr = dxdr(p0, p1, p2, p3, p4, p5, p6, p7, &[r0, r1, r2]);
            use slice_of_array::SliceFlatExt;
            let detjac =
                crate::mat3_col_major::determinant::<Real>(dxdr.flat().try_into().unwrap());
            detjac * w
        })
        .sum();
    vol
}

#[test]
fn test_volume() {
    let node2xyz: [[f64; 3]; 8] = [
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [1., 0., 1.],
        [1., 1., 1.],
        [0., 1., 1.],
    ];
    let v = volume(
        &node2xyz[0],
        &node2xyz[1],
        &node2xyz[2],
        &node2xyz[3],
        &node2xyz[4],
        &node2xyz[5],
        &node2xyz[6],
        &node2xyz[7],
        2,
    );
    assert!((v - 1.0).abs() < 1.0e-10, "volume of unit cube: {}", v);
}

/// Returns parametric coordinates (r0,r1,r2) of the origin inside the hex, or None if outside.
/// Valid domain: r0,r1,r2 in [-1,1].
/// Vertices: [xyz, Xyz, XYz, xYz, xyZ, XyZ, XYZ, xYZ] (lowercase=-1, uppercase=+1).
pub fn parametric_coord_for_origin<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    p6: &[Real; 3],
    p7: &[Real; 3],
) -> Option<[Real; 3]>
where
    Real: num_traits::Float,
{
    use slice_of_array::SliceFlatExt;
    let zero = Real::zero();
    let one = Real::one();
    let mut pco = [zero; 3];
    for _ in 0..20 {
        let an = shapefunc(&pco);
        let pos = [
            an[0] * p0[0]
                + an[1] * p1[0]
                + an[2] * p2[0]
                + an[3] * p3[0]
                + an[4] * p4[0]
                + an[5] * p5[0]
                + an[6] * p6[0]
                + an[7] * p7[0],
            an[0] * p0[1]
                + an[1] * p1[1]
                + an[2] * p2[1]
                + an[3] * p3[1]
                + an[4] * p4[1]
                + an[5] * p5[1]
                + an[6] * p6[1]
                + an[7] * p7[1],
            an[0] * p0[2]
                + an[1] * p1[2]
                + an[2] * p2[2]
                + an[3] * p3[2]
                + an[4] * p4[2]
                + an[5] * p5[2]
                + an[6] * p6[2]
                + an[7] * p7[2],
        ];
        let jac = dxdr(p0, p1, p2, p3, p4, p5, p6, p7, &pco);
        let jac_col: &[Real; 9] = jac.flat().try_into().unwrap();
        let jac_col = crate::mat3_col_major::transpose(jac_col);
        let j_inv = crate::mat3_col_major::try_inverse(&jac_col)?;
        let dpco = crate::mat3_col_major::mult_vec(&j_inv, &pos);
        for i in 0..3 {
            pco[i] = pco[i] - dpco[i];
        }
    }
    if pco.iter().all(|&x| x >= -one && x <= one) {
        Some(pco)
    } else {
        None
    }
}

/// Returns the nearest point to the origin on the hex and its parametric coordinates (r0,r1,r2).
/// Valid domain: r0,r1,r2 in [-1,1].
/// Vertices: [xyz, Xyz, XYz, xYz, xyZ, XyZ, XYZ, xYZ] (lowercase=-1, uppercase=+1).
pub fn nearest_to_origin<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    p3: &[Real; 3],
    p4: &[Real; 3],
    p5: &[Real; 3],
    p6: &[Real; 3],
    p7: &[Real; 3],
) -> ([Real; 3], [Real; 3])
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;

    if let Some(pco) = parametric_coord_for_origin(p0, p1, p2, p3, p4, p5, p6, p7) {
        return ([zero; 3], pco);
    }

    // p0 is at (r0=-1, r1=-1, r2=-1)
    let mut p_min = *p0;
    let mut pco_min = [-one, -one, -one];
    let mut d_min = crate::vec3::norm(p0);

    let mut update = |p: [Real; 3], pco: [Real; 3]| {
        let d = crate::vec3::norm(&p);
        if d < d_min {
            d_min = d;
            p_min = p;
            pco_min = pco;
        }
    };

    // 6 quad faces: quad3 returns (s0,s1) in [0,1]^2; hex uses [-1,1]^3, so r = 2s-1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p0, p1, p2, p3);
        update(p, [two * s0 - one, two * s1 - one, -one]);
    } // r2=-1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p4, p5, p6, p7);
        update(p, [two * s0 - one, two * s1 - one, one]);
    } // r2=+1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p0, p1, p5, p4);
        update(p, [two * s0 - one, -one, two * s1 - one]);
    } // r1=-1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p3, p2, p6, p7);
        update(p, [two * s0 - one, one, two * s1 - one]);
    } // r1=+1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p0, p3, p7, p4);
        update(p, [-one, two * s0 - one, two * s1 - one]);
    } // r0=-1
    {
        let (p, s0, s1) = crate::quad3::nearest_to_origin(p1, p2, p6, p5);
        update(p, [one, two * s0 - one, two * s1 - one]);
    } // r0=+1

    // 12 edges: edge3 returns (p, s0=1-t, s1=t); r of free axis = 2*s1-1 or 2*s0-1 depending on direction
    // bottom face (r2=-1)
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p0, p1);
        update(p, [two * s1 - one, -one, -one]);
    } // r0 free
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p1, p2);
        update(p, [one, two * s1 - one, -one]);
    } // r1 free
    {
        let (p, s0, _s1) = crate::edge3::nearest_to_origin3(p2, p3);
        update(p, [two * s0 - one, one, -one]);
    } // r0 free (reversed)
    {
        let (p, s0, _s1) = crate::edge3::nearest_to_origin3(p3, p0);
        update(p, [-one, two * s0 - one, -one]);
    } // r1 free (reversed)
    // top face (r2=+1)
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p4, p5);
        update(p, [two * s1 - one, -one, one]);
    }
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p5, p6);
        update(p, [one, two * s1 - one, one]);
    }
    {
        let (p, s0, _s1) = crate::edge3::nearest_to_origin3(p6, p7);
        update(p, [two * s0 - one, one, one]);
    }
    {
        let (p, s0, _s1) = crate::edge3::nearest_to_origin3(p7, p4);
        update(p, [-one, two * s0 - one, one]);
    }
    // vertical edges (r2 free)
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p0, p4);
        update(p, [-one, -one, two * s1 - one]);
    }
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p1, p5);
        update(p, [one, -one, two * s1 - one]);
    }
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p2, p6);
        update(p, [one, one, two * s1 - one]);
    }
    {
        let (p, _s0, s1) = crate::edge3::nearest_to_origin3(p3, p7);
        update(p, [-one, one, two * s1 - one]);
    }

    (p_min, pco_min)
}

#[test]
fn test_nearest_to_origin() {
    // origin inside: hex from [-2,-2,-2] to [2,2,2]
    {
        let p0 = [-2.0f64, -2.0, -2.0];
        let p1 = [2.0, -2.0, -2.0];
        let p2 = [2.0, 2.0, -2.0];
        let p3 = [-2.0, 2.0, -2.0];
        let p4 = [-2.0, -2.0, 2.0];
        let p5 = [2.0, -2.0, 2.0];
        let p6 = [2.0, 2.0, 2.0];
        let p7 = [-2.0, 2.0, 2.0];
        let (p, pco) = nearest_to_origin(&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7);
        assert!(
            crate::vec3::norm(&p) < 1.0e-10,
            "inside: expected origin, got {:?}",
            p
        );
        let (r0, r1, r2) = (pco[0], pco[1], pco[2]);
        assert!(
            r0 > -1.0 && r0 < 1.0 && r1 > -1.0 && r1 < 1.0 && r2 > -1.0 && r2 < 1.0,
            "inside: pco should be interior: {:?}",
            pco
        );
        let an = shapefunc(&pco);
        let verts = [p0, p1, p2, p3, p4, p5, p6, p7];
        let recon: [f64; 3] =
            std::array::from_fn(|i| an.iter().zip(verts.iter()).map(|(&w, v)| w * v[i]).sum());
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction failed"
        );
    }
    // origin outside: hex shifted to x in [3,7]
    {
        let p0 = [3.0f64, -1.0, -1.0];
        let p1 = [7.0, -1.0, -1.0];
        let p2 = [7.0, 1.0, -1.0];
        let p3 = [3.0, 1.0, -1.0];
        let p4 = [3.0, -1.0, 1.0];
        let p5 = [7.0, -1.0, 1.0];
        let p6 = [7.0, 1.0, 1.0];
        let p7 = [3.0, 1.0, 1.0];
        let (p, pco) = nearest_to_origin(&p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7);
        let (r0, r1, r2) = (pco[0], pco[1], pco[2]);
        assert!(
            r0 >= -1.0 - 1.0e-10
                && r0 <= 1.0 + 1.0e-10
                && r1 >= -1.0 - 1.0e-10
                && r1 <= 1.0 + 1.0e-10
                && r2 >= -1.0 - 1.0e-10
                && r2 <= 1.0 + 1.0e-10,
            "outside: pco out of domain: {:?}",
            pco
        );
        let an = shapefunc(&pco);
        let verts = [p0, p1, p2, p3, p4, p5, p6, p7];
        let recon: [f64; 3] =
            std::array::from_fn(|i| an.iter().zip(verts.iter()).map(|(&w, v)| w * v[i]).sum());
        assert!(
            crate::vec3::norm(&crate::vec3::sub(&recon, &p)) < 1.0e-10,
            "reconstruction failed"
        );
        let d = crate::vec3::norm(&p);
        for v in &verts {
            assert!(
                d <= crate::vec3::norm(v) + 1.0e-10,
                "vertex {:?} is closer than result",
                v
            );
        }
    }
}

pub fn iso_surface(
    tri2xyz: &mut Vec<[f64; 3]>,
    center: [f64; 3],
    half_width: f64,
    corner_value: &[f64; 8],
) {
    // Evaluate Dis Node
    let zero = 0f64;
    let cubeindex = {
        let mut cubeindex = 0;
        for i_node_hex in 0..8 {
            if corner_value[i_node_hex] < zero {
                cubeindex |= 1 << i_node_hex;
            };
        }
        cubeindex
    };
    // dbg!(cubeindex, cent_);

    let vert_list = {
        let mut vertlist = [[f64::INFINITY; 3]; 12];
        let e0 = ISO_EDGE[cubeindex];
        for i_edge in 0..12 {
            if (e0 >> i_edge) & 1 == 0 {
                continue;
            }
            let i0_node = EDGE2NODE[i_edge][0];
            let i1_node = EDGE2NODE[i_edge][1];
            iso_position_on_edge(
                &mut vertlist[i_edge],
                center,
                half_width,
                i0_node,
                i1_node,
                corner_value[i0_node],
                corner_value[i1_node],
            );
            // println!("  {} {} {:?} {:?} {}", i0_node, i1_node, &vertlist[i_edge], cent_, hw_);
        }
        vertlist
    };

    let ti = &ISO_TRI[cubeindex];
    for tri in ti.chunks(3) {
        if tri[0] == -1 {
            break;
        }
        assert_eq!(tri.len(), 3);
        tri2xyz.push(vert_list[tri[0] as usize]);
        tri2xyz.push(vert_list[tri[1] as usize]);
        tri2xyz.push(vert_list[tri[2] as usize]);
    }
}

pub fn subdivide<INDEX>(
    corner: &[INDEX; 8],
    edge: &[INDEX; 12],
    quad: &[INDEX; 6],
    center: INDEX,
) -> [[INDEX; 8]; 8]
where
    INDEX: num_traits::PrimInt,
{
    let (c, e, q, ct) = (corner, edge, quad, center);
    // Each sub-hex i has corner[i] at its local "origin" corner, growing toward ct.
    // Node ordering: [xyz, Xyz, XYz, xYz, xyZ, XyZ, XYZ, xYZ]
    // quad: [0]=r2=-1(bottom), [1]=r1=-1(front), [2]=r0=+1(right),
    //        [3]=r1=+1(back),  [4]=r0=-1(left),  [5]=r2=+1(top)
    // edge: [0]=0-1, [1]=1-2, [2]=2-3, [3]=3-0, [4]=4-5, [5]=5-6,
    //        [6]=6-7, [7]=7-4, [8]=0-4, [9]=1-5, [10]=2-6, [11]=3-7
    [
        [c[0], e[0], q[0], e[3], e[8], q[1], ct, q[4]], // (−1,−1,−1)
        [e[0], c[1], e[1], q[0], q[1], e[9], q[2], ct], // (+1,−1,−1)
        [q[0], e[1], c[2], e[2], ct, q[2], e[10], q[3]], // (+1,+1,−1)
        [e[3], q[0], e[2], c[3], q[4], ct, q[3], e[11]], // (−1,+1,−1)
        [e[8], q[1], ct, q[4], c[4], e[4], q[5], e[7]], // (−1,−1,+1)
        [q[1], e[9], q[2], ct, e[4], c[5], e[5], q[5]], // (+1,−1,+1)
        [ct, q[2], e[10], q[3], q[5], e[5], c[6], e[6]], // (+1,+1,+1)
        [q[4], ct, q[3], e[11], e[7], q[5], e[6], c[7]], // (−1,+1,+1)
    ]
}
