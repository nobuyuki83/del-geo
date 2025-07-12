#![allow(clippy::needless_range_loop)]
/*
vec < mat < aabb < obb <
line < ray < edge <
plane < tri < quad <
sphere < tet
*/

pub mod aabb;
pub mod aabb2;
pub mod aabb3;
pub mod range;

pub mod curve_linear_coords;
pub mod mat2x3_col_major;

pub mod mat3_array_of_array;
pub mod mat3_col_major;
pub mod mat4_col_major;
pub mod obb3;
pub mod vec2;
pub mod vec3;
//
pub mod bezier_cubic;
pub mod bezier_quadratic;
pub mod ccd2;
pub mod ccd3;
pub mod edge;
pub mod edge2;
pub mod edge3;
pub mod hex;
pub mod line2;
pub mod mat2_col_major;
pub mod mat2_sym;
pub mod mat3_row_major;
pub mod mat3_sym;
pub mod mat3x4_col_major;
pub mod matn_row_major;
pub mod obb2;
pub mod plane;
pub mod polynomial_root;
pub mod quaternion;
pub mod sphere;
pub mod spherical_harmonics;
pub mod tet;
pub mod tri2;
pub mod tri3;
pub mod uvec3;
pub mod vecn;
pub mod view_projection;
pub mod view_rotation;
pub mod ndc;
