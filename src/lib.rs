
/*
 aabb ->
 line -> ray -> edge ->
 polyline -> curve_quadratic -> curve_cubic -> curve_ndegree ->
 plane -> tri -> quad ->
 tet
 */
pub mod vec3;
//
pub mod line;
pub mod line3;
pub mod edge;
pub mod edge2;
pub mod edge3;
//
pub mod polyline;
//
pub mod plane;
pub mod tri2;
pub mod tri3;
pub mod polygon;
//
pub mod tet;

