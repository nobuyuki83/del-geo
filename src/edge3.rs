use num_traits::AsPrimitive;

pub fn nearest_point3_<T>(
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
        if vec3::dot_(&d, &d) > 1.0e-20_f64.as_() {
            let ps = [
                edge_pos0[0] - point_pos[0],
                edge_pos0[1] - point_pos[1],
                edge_pos0[2] - point_pos[2]];
            let a = vec3::dot_(&d, &d);
            let b = vec3::dot_(&d, &ps);
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

////////////////////////////////////

pub fn nearest_to_line3(
    edge_start: &nalgebra::Vector3<f64>,
    edge_end: &nalgebra::Vector3<f64>,
    line_origin: &nalgebra::Vector3<f64>,
    line_direction: &nalgebra::Vector3<f64>) -> (nalgebra::Vector3<f64>, nalgebra::Vector3<f64>)
{
    let (scale, scaled_ratio_edge, _, scaled_nearest_edge, scaled_nearest_line)
        = crate::line3::nearest_to_line3(
        edge_start, &(edge_end - edge_start),
        line_origin, line_direction);
    if scale.abs() < 1.0e-10 { // pararell
        let nearest_edge = (edge_start + edge_end) * 0.5;
        let (nearest_line, _) = crate::line::nearest_to_point(
            &nearest_edge, &line_origin, &line_direction);
        return (nearest_edge, nearest_line);
    }
    let ratio_edge = scaled_ratio_edge / scale;
    if ratio_edge > 0_f64 && ratio_edge < 1_f64 { // nearst point is inside the segment
        let nearest_edge = scaled_nearest_edge / scale;
        let nearest_line = scaled_nearest_line / scale;
        return (nearest_edge, nearest_line);
    }
    //
    let (p1, _) = crate::line::nearest_to_point(edge_start, line_origin, line_direction);
    let (p2, _) = crate::line::nearest_to_point(edge_end, line_origin, line_direction);
    let dist1 = (p1 - edge_start).norm();
    let dist2 = (p2 - edge_end).norm();
    if dist1 < dist2 {
        let nearest_edge = edge_start;
        let nearest_line = p1;
        return (*nearest_edge, nearest_line);
    }
    let nearest_edge = edge_end;
    let nearest_line = p2;
    return (*nearest_edge, nearest_line);
}


