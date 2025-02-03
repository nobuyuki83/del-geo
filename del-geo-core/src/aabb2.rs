//! methods for 2D Axis-aligned Bounding Box (AABB)

pub fn from_point<T>(p: &[T; 2], rad: T) -> [T; 4]
where
    T: num_traits::Float,
{
    assert!(rad >= T::zero());
    [p[0] - rad, p[1] - rad, p[0] + rad, p[1] + rad]
}

pub fn from_two_points<T>(p0: &[T; 2], p1: &[T; 2], rad: T) -> [T; 4]
where
    T: num_traits::Float,
{
    [
        (p0[0] - rad).min(p1[0] - rad),
        (p0[1] - rad).min(p1[1] - rad),
        (p0[0] + rad).max(p1[0] + rad),
        (p0[1] + rad).max(p1[1] + rad),
    ]
}

pub fn from_two_aabbs<T>(i0: &[T; 4], i1: &[T; 4]) -> [T; 4]
where
    T: num_traits::Float,
{
    let mut o = [T::zero(); 4];
    for i in 0..2 {
        o[i] = i0[i].min(i1[i]);
        o[i + 2] = i0[i + 2].max(i1[i + 2]);
    }
    o
}

// above: from method
// -----------------------

pub fn add_point<T>(aabb2: &mut [T; 4], &[p0, p1]: &[T; 2], rad: T)
where
    T: num_traits::Float,
{
    assert!(rad >= T::zero());
    aabb2[0] = aabb2[0].min(p0 - rad);
    aabb2[1] = aabb2[1].min(p1 - rad);
    aabb2[2] = aabb2[2].max(p0 + rad);
    aabb2[3] = aabb2[3].max(p1 + rad);
}

pub fn rasterize<T>(aabb: &[T; 4], img_size: &(usize, usize)) -> [usize; 4]
where
    T: num_traits::Float + num_traits::AsPrimitive<usize> + 'static + Copy,
    usize: num_traits::AsPrimitive<T>,
{
    use num_traits::AsPrimitive;
    let half = T::one() / (T::one() + T::one());
    [
        (aabb[0] - half).ceil().max(T::zero()).as_(),
        (aabb[1] - half).ceil().max(T::zero()).as_(),
        (aabb[2] - half).ceil().min(img_size.0.as_()).as_(),
        (aabb[3] - half).ceil().min(img_size.1.as_()).as_(),
    ]
}

pub fn center<T>(aabb: &[T; 4]) -> [T; 2]
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    [(aabb[0] + aabb[2]) * half, (aabb[1] + aabb[3]) * half]
}

pub fn max_edge_size<T>(aabb: &[T; 4]) -> T
where
    T: num_traits::Float,
{
    let lx = aabb[2] - aabb[0];
    let ly = aabb[3] - aabb[1];
    lx.max(ly)
}

pub fn transform_homogeneous<T>(aabb: &[T; 4], mat3_col_major: &[T; 9]) -> [T; 4]
where
    T: num_traits::Float,
{
    use crate::mat3_col_major::Mat3ColMajor;
    let p0 = mat3_col_major
        .transform_homogeneous(&[aabb[0], aabb[1]])
        .unwrap();
    let p1 = mat3_col_major
        .transform_homogeneous(&[aabb[0], aabb[3]])
        .unwrap();
    let p2 = mat3_col_major
        .transform_homogeneous(&[aabb[2], aabb[1]])
        .unwrap();
    let p3 = mat3_col_major
        .transform_homogeneous(&[aabb[2], aabb[3]])
        .unwrap();
    let ax = [p0[0], p1[0], p2[0], p3[0]];
    let ay = [p0[1], p1[1], p2[1], p3[1]];
    [
        ax.iter().fold(T::infinity(), |a, &b| a.min(b)),
        ay.iter().fold(T::infinity(), |a, &b| a.min(b)),
        ax.iter().fold(T::neg_infinity(), |a, &b| a.max(b)),
        ay.iter().fold(T::neg_infinity(), |a, &b| a.max(b)),
    ]
}

pub fn sample<Reng, T>(aabb: &[T; 4], reng: &mut Reng) -> [T; 2]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<T>,
{
    let r0 = reng.random::<T>();
    let r1 = reng.random::<T>();
    [
        aabb[0] + r0 * (aabb[2] - aabb[0]),
        aabb[1] + r1 * (aabb[3] - aabb[1]),
    ]
}

/// transform aabb to unit square (0,1)^2 while preserving aspect ratio
/// return 3x3 homogeneous transformation matrix in **column major** order
pub fn to_transformation_world2unit_ortho_preserve_asp(aabb_world: &[f32; 4]) -> [f32; 9] {
    let cntr = center(aabb_world);
    let size = max_edge_size(aabb_world);
    let a = 1f32 / size;
    let b = -a * cntr[0] + 0.5;
    let c = -a * cntr[1] + 0.5;
    [a, 0., 0., 0., a, 0., b, c, 1.]
}

#[test]
fn test_to_transformation_world2unit_ortho_preserve_asp() {
    let aabb = [-4.1, 0.3, 4.5, 3.3];
    let transf = to_transformation_world2unit_ortho_preserve_asp(&aabb);
    let mut reng = rand::rng();
    for _ in 0..1000 {
        let p0 = sample(&aabb, &mut reng);
        let q0 = crate::mat3_col_major::transform_homogeneous(&transf, &p0).unwrap();
        assert!(q0[0] >= 0.0 && q0[0] <= 1.0);
        assert!(q0[1] >= 0.0 && q0[1] <= 1.0);
    }
}

pub fn is_include_point2<Real>(aabb: &[Real; 4], &[p0, p1]: &[Real; 2]) -> bool
where
    Real: num_traits::Float,
{
    p0 >= aabb[0] && p0 <= aabb[2] && p1 >= aabb[1] && p1 <= aabb[3]
}

/// signed distance from axis-aligned bounding box
/// * `pos_in` - where the signed distance is evaluated
/// * `x_min` - bounding box's x-coordinate minimum
/// * `x_max` - bounding box's x-coordinate maximum
/// * `y_min` - bounding box's y-coordinate minimum
/// * `y_max` - bounding box's y-coordinate maximum
/// * signed distance (inside is negative)
pub fn is_intersect_square<Real>(aabb: &[Real; 4], pos_in: &[Real; 2], rad: Real) -> bool
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let half = one / (one + one);
    let min0 = [aabb[0], aabb[1]];
    let max0 = [aabb[2], aabb[3]];
    let x_center = (max0[0] + min0[0]) * half;
    let y_center = (max0[1] + min0[1]) * half;
    let x_dist = (pos_in[0] - x_center).abs() - (max0[0] - min0[0]) * half;
    let y_dist = (pos_in[1] - y_center).abs() - (max0[1] - min0[1]) * half;
    rad > x_dist.max(y_dist)
}

pub fn nearest_point2(aabb: &[f32; 4], &[p0, p1]: &[f32; 2]) -> [f32; 2] {
    [p0.clamp(aabb[0], aabb[2]), p1.clamp(aabb[1], aabb[3])]
}

/// <https://iquilezles.org/articles/distfunctions2d/>
pub fn sdf<Real>(aabb: &[Real; 4], p: &[Real; 2]) -> Real
where
    Real: num_traits::Float,
{
    let zero = Real::zero();
    let one = Real::one();
    let half = one / (one + one);
    let c = center(aabb);
    let b = [(aabb[2] - aabb[0]) * half, (aabb[3] - aabb[1]) * half];
    let p = [p[0] - c[0], p[1] - c[1]];
    let d = [p[0].abs() - b[0], p[1].abs() - b[1]];
    crate::vec2::length(&[d[0].max(zero), d[1].max(zero)]) + d[0].max(d[1]).min(zero)
}

/*
/// signed distance from axis-aligned bounding box
/// * `pos_in` - where the signed distance is evaluated
/// * `x_min` - bounding box's x-coordinate minimum
/// * `x_max` - bounding box's x-coordinate maximum
/// * `y_min` - bounding box's y-coordinate minimum
/// * `y_max` - bounding box's y-coordinate maximum
/// * signed distance (inside is negative)
pub fn signed_distance<Real>(
    pos_in: nalgebra::Vector2<Real>,
    min0: nalgebra::Vector2<Real>,
    max0: nalgebra::Vector2<Real>,
) -> Real
where
    Real: nalgebra::RealField + Copy,
    f64: AsPrimitive<Real>,
{
    let half = 0.5_f64.as_();
    let x_center = (max0.x + min0.x) * half;
    let y_center = (max0.y + min0.y) * half;
    let x_dist = (pos_in.x - x_center).abs() - (max0.x - min0.x) * half;
    let y_dist = (pos_in.y - y_center).abs() - (max0.y - min0.y) * half;
    x_dist.max(y_dist)
}
 */

#[test]
pub fn test_sdf() {
    let aabb: [f64; 4] = [0.2, 0.4, 3.0, 4.0];
    assert!((sdf(&aabb, &[0.2, 0.4])) < f32::EPSILON as f64);
    assert!((sdf(&aabb, &[0.0, 0.0]) - (0.2f64 * 0.2 + 0.4 * 0.4).sqrt()).abs() < f64::EPSILON);
    assert!((sdf(&aabb, &[1.6, 2.2]) + 1.4) < f64::EPSILON);
    assert!(sdf(&aabb, &[0.2, 2.0]).abs() < f32::EPSILON as f64);
    assert!((sdf(&aabb, &[0.4, 2.0]) + 0.2).abs() < f32::EPSILON as f64);
    assert!((sdf(&aabb, &[0.0, 2.0]) - 0.2).abs() < f32::EPSILON as f64);
}

pub fn scale<Real>(aabb: &[Real; 4], s: Real) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let one = Real::one();
    let half = one / (one + one);
    [
        (aabb[0] + aabb[2]) * half - (aabb[2] - aabb[0]) * half * s,
        (aabb[1] + aabb[3]) * half - (aabb[3] - aabb[1]) * half * s,
        (aabb[0] + aabb[2]) * half + (aabb[2] - aabb[0]) * half * s,
        (aabb[1] + aabb[3]) * half + (aabb[3] - aabb[1]) * half * s,
    ]
}

pub fn translate<Real>(aabb: &[Real; 4], t: &[Real; 2]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    [
        aabb[0] + t[0],
        aabb[1] + t[1],
        aabb[2] + t[0],
        aabb[3] + t[1],
    ]
}

pub fn overlapping_tiles(
    aabb: &[f32; 4],
    tile_size: usize,
    tile_shape: (usize, usize),
) -> std::collections::BTreeSet<usize> {
    let ix0 = (aabb[0] / tile_size as f32).floor() as i32;
    let iy0 = (aabb[1] / tile_size as f32).floor() as i32;
    let ix1 = (aabb[2] / tile_size as f32).floor() as i32 + 1;
    let iy1 = (aabb[3] / tile_size as f32).floor() as i32 + 1;
    let mut tiles = std::collections::BTreeSet::<usize>::new();
    for ix in ix0..ix1 {
        assert_ne!(ix, ix1);
        if ix < 0 || ix >= (tile_shape.0 as i32) {
            continue;
        }
        let ix = ix as usize;
        for iy in iy0..iy1 {
            assert_ne!(iy, iy1);
            if iy < 0 || iy >= (tile_shape.1 as i32) {
                continue;
            }
            let iy = iy as usize;
            let i_tile = iy * tile_shape.0 + ix;
            tiles.insert(i_tile);
        }
    }
    tiles
}

// -------------------------------------------------------------------------------

pub type AABB2<Real> = crate::aabb::AABB<Real, 2, 4>;

pub fn from_slice<Real>(s: &[Real]) -> AABB2<Real>
where
    Real: num_traits::Float,
{
    AABB2 {
        aabb: s[..4].try_into().unwrap(),
    }
}

pub fn from_aabbs<Real>(aabbs: &[Real], i_aabb: usize) -> AABB2<Real>
where
    Real: num_traits::Float,
{
    AABB2 {
        aabb: aabbs[(i_aabb * 4)..(i_aabb * 4) + 4].try_into().unwrap(),
    }
}
