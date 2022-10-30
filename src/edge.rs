pub fn nearest_edge3_point3(
    point_pos: &[f32],
    edge_pos0: &[f32],
    edge_pos1: &[f32]) -> [f32; 3]
{
    use crate::vec3;
    let d = [
        edge_pos1[0] - edge_pos0[0],
        edge_pos1[1] - edge_pos0[1],
        edge_pos1[2] - edge_pos0[2]];
    let t = {
        if vec3::dot(&d, &d) > 1.0e-20 {
            let ps = [
                edge_pos0[0] - point_pos[0],
                edge_pos0[1] - point_pos[1],
                edge_pos0[2] - point_pos[2]];
            let a = vec3::dot(&d, &d);
            let b = vec3::dot(&d, &ps);
            (-b / a).clamp(0_f32, 1_f32)
        } else {
            0.5_f32
        }
    };
    [edge_pos0[0] + t * d[0],
        edge_pos0[1] + t * d[1],
        edge_pos0[2] + t * d[2]]
}