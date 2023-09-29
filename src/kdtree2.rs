use num_traits::Zero;

#[derive(Clone)]
pub struct Node<Real> {
    pos: (nalgebra::Vector2::<Real>, usize),
    idx_node_left: usize,
    idx_node_right: usize,
}

impl<Real> Node<Real>
    where Real: num_traits::Zero
{
    pub fn new() -> Self {
        Node {
            pos: (nalgebra::Vector2::<Real>::new(Real::zero(), Real::zero()), usize::MAX),
            idx_node_left: usize::MAX,
            idx_node_right: usize::MAX,
        }
    }
}


/// construct Kd-tree recursively
/// * `nodes`
/// * `idx_node`
/// * `points`
/// * `idx_point_begin`
/// * `idx_point_end`
/// * `i_depth`
pub fn construct_kdtree<Real>(
    nodes: &mut Vec<Node<Real>>,
    idx_node: usize,
    points: &mut Vec<(nalgebra::Vector2::<Real>, usize)>,
    idx_point_begin: usize,
    idx_point_end: usize,
    i_depth: i32)
    where Real: nalgebra::RealField + Copy
{
    if idx_point_end - idx_point_begin == 1 { // leaf node
        nodes[idx_node].pos = points[idx_point_begin];
        return;
    }

    if i_depth % 2 == 0 { // sort by x-coordinate
        points[idx_point_begin..idx_point_end].sort_by(
            |a, b| a.0.x.partial_cmp(&b.0.x).unwrap());
    } else { // sort by y-coordinate
        points[idx_point_begin..idx_point_end].sort_by(
            |a, b| a.0.y.partial_cmp(&b.0.y).unwrap());
    }

    let idx_point_mid = (idx_point_end - idx_point_begin) / 2 + idx_point_begin; // median point
    nodes[idx_node].pos = points[idx_point_mid];

    if idx_point_begin != idx_point_mid { // there is at least one point smaller than median
        let idx_node_left = nodes.len();
        nodes.resize(nodes.len() + 1, Node::new());
        nodes[idx_node].idx_node_left = idx_node_left;
        construct_kdtree(
            nodes, idx_node_left,
            points, idx_point_begin, idx_point_mid,
            i_depth + 1);
    }
    if idx_point_mid + 1 != idx_point_end { // there is at least one point larger than median
        let idx_node_right = nodes.len();
        nodes.resize(nodes.len() + 1, Node::new());
        nodes[idx_node].idx_node_right = idx_node_right;
        construct_kdtree(
            nodes, idx_node_right,
            points, idx_point_mid + 1, idx_point_end,
            i_depth + 1);
    }
}

fn find_edges<Real>(
    xyz: &mut Vec<Real>,
    nodes: &Vec<Node<Real>>,
    idx_node: usize,
    min: nalgebra::Vector2::<Real>,
    max: nalgebra::Vector2::<Real>,
    i_depth: i32)
where Real: Copy
{
    if idx_node >= nodes.len() { return; }
    let pos = &nodes[idx_node].pos.0;
    if i_depth % 2 == 0 {
        xyz.push(pos[0]);
        xyz.push(min[1]);
        xyz.push(pos[0]);
        xyz.push(max[1]);
        find_edges(xyz, nodes, nodes[idx_node].idx_node_left,
                          min,
                          nalgebra::Vector2::new(pos[0], max[1]),
                          i_depth + 1);
        find_edges(xyz, nodes, nodes[idx_node].idx_node_right,
                          nalgebra::Vector2::new(pos[0], min[1]),
                          max,
                          i_depth + 1);
    } else {
        xyz.push(min[0]);
        xyz.push(pos[1]);
        xyz.push(max[0]);
        xyz.push(pos[1]);
        find_edges(xyz, nodes, nodes[idx_node].idx_node_left,
                          min,
                          nalgebra::Vector2::new(max[0], pos[1]),
                          i_depth + 1);
        find_edges(xyz, nodes, nodes[idx_node].idx_node_right,
                          nalgebra::Vector2::new(min[0], pos[1]),
                          max,
                          i_depth + 1);
    }
}


pub struct KdTree2<Real> {
    min: nalgebra::Vector2::<Real>,
    max: nalgebra::Vector2::<Real>,
    nodes: Vec<Node<Real>>,
}

impl<Real> KdTree2<Real>
    where Real: num_traits::Float + nalgebra::RealField
{
    pub fn from_matrix(points_: &nalgebra::Matrix2xX::<Real>) -> Self {
        let mut ps = points_
            .column_iter().enumerate()
            .map(|v| (v.1.into_owned(), v.0))
            .collect();
        let mut nodes = Vec::<Node<Real>>::new();
        nodes.resize(1, Node::new());
        construct_kdtree(&mut nodes, 0,
                         &mut ps, 0, points_.shape().1,
                         0);
        use num_traits::Float;
        let min_x = points_.column_iter().fold(Real::nan(), |m, v| Float::min(v.x, m));
        let max_x = points_.column_iter().fold(Real::nan(), |m, v| Float::max(v.x, m));
        let min_y = points_.column_iter().fold(Real::nan(), |m, v| Float::min(v.y, m));
        let max_y = points_.column_iter().fold(Real::nan(), |m, v| Float::max(v.y, m));
        KdTree2 {
            min: nalgebra::Vector2::<Real>::new(min_x, min_y),
            max: nalgebra::Vector2::<Real>::new(max_x, max_y),
            nodes: nodes,
        }
    }

    pub fn edges(&self) -> Vec<Real> {
        let mut xys = Vec::<Real>::new();
        find_edges(
            &mut xys,
            &self.nodes,
            0,
            self.min, self.max,
            0);
        xys.push(self.min.x);
        xys.push(self.min.y);
        xys.push(self.max.x);
        xys.push(self.min.y);
        //
        xys.push(self.max.x);
        xys.push(self.min.y);
        xys.push(self.max.x);
        xys.push(self.max.y);
        //
        xys.push(self.max.x);
        xys.push(self.max.y);
        xys.push(self.min.x);
        xys.push(self.max.y);
        //
        xys.push(self.min.x);
        xys.push(self.max.y);
        xys.push(self.min.x);
        xys.push(self.min.y);
        xys
    }
}