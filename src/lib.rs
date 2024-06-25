/*
vec -> mat -> aabb ->
line -> ray -> edge ->
plane -> tri -> quad ->
tet
*/

pub mod aabb;
pub mod aabb2;
pub mod aabb3;
pub mod mat2;
pub mod mat3;
pub mod mat4;
pub mod vec2;
pub mod vec3;
//
pub mod bezier_cubic;
pub mod bezier_quadratic;
pub mod ccd;
pub mod edge;
pub mod edge2;
pub mod edge3;
pub mod line;
pub mod line2;
pub mod line3;
mod matn;
pub mod plane;
pub mod sphere;
pub mod tet;
pub mod tri2;
pub mod tri3;
