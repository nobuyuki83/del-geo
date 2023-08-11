use num_traits::AsPrimitive;

pub fn nearest_point3<T>(
    point_pos: &[T],
    edge_pos0: &[T],
    edge_pos1: &[T]) -> [T; 3]
    where T: num_traits::Float + 'static + Copy + PartialOrd,
          f64: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    let d = [
        edge_pos1[0] - edge_pos0[0],
        edge_pos1[1] - edge_pos0[1],
        edge_pos1[2] - edge_pos0[2]];
    let t = {
        if vec3::dot(&d, &d) > 1.0e-20_f64.as_() {
            let ps = [
                edge_pos0[0] - point_pos[0],
                edge_pos0[1] - point_pos[1],
                edge_pos0[2] - point_pos[2]];
            let a = vec3::dot(&d, &d);
            let b = vec3::dot(&d, &ps);
            let mut r: T = -b / a;
            r = if r < 0_f64.as_() { 0_f64.as_() } else { r };
            r = if r > 1_f64.as_() { 1_f64.as_() } else { r };
            r
        } else {
            0.5_f64.as_()
        }
    };
    [edge_pos0[0] + t * d[0],
        edge_pos0[1] + t * d[1],
        edge_pos0[2] + t * d[2]]
}


