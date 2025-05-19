type Vec2f = nalgebra::Vector2<f32>;

/// node of the Kd-tree
/// Kd-tree is a binary tree. Each node has left and right children nodes
#[derive(Debug, Clone)]
pub struct Node {
    pub pos: Vec2f,
    pub idx_node_left: usize,
    pub idx_node_right: usize,
}

impl std::default::Default for Node {
    fn default() -> Self {
        Self {
            pos: Vec2f::new(0., 0.),
            idx_node_left: usize::MAX,
            idx_node_right: usize::MAX,
        }
    }
}

impl Node {
    pub fn new(x: f32, y: f32) -> Self {
        Self {
            pos: Vec2f::new(x, y),
            idx_node_left: usize::MAX,
            idx_node_right: usize::MAX,
        }
    }
}

/// construction of the Kd-tree
pub fn construct_kdtree(
    nodes: &mut Vec<Node>,
    idx_node: usize,
    points: &mut [Vec2f],
    idx_point_begin: usize,
    idx_point_end: usize,
    i_depth: usize,
) {
    if idx_point_end - idx_point_begin == 1 {
        nodes[idx_node].pos = points[idx_point_begin];
        return;
    }

    if i_depth % 2 == 0 {
        // if depth is even, split in the horizontal direction
        points[idx_point_begin..idx_point_end].sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());
    } else {
        // if depth is even, split in the vertical direction
        points[idx_point_begin..idx_point_end].sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
    }

    // split in the middle
    let idx_point_mid = (idx_point_end - idx_point_begin) / 2 + idx_point_begin;
    nodes[idx_node].pos = points[idx_point_mid];

    if idx_point_begin != idx_point_mid {
        let idx_node_left = nodes.len();
        nodes.push(Node::new(0.0, 0.0));
        nodes[idx_node].idx_node_left = idx_node_left;
        construct_kdtree(
            nodes,
            idx_node_left,
            points,
            idx_point_begin,
            idx_point_mid,
            i_depth + 1,
        );
    }

    if idx_point_mid + 1 != idx_point_end {
        let idx_node_right = nodes.len();
        nodes.push(Node::new(0.0, 0.0));
        nodes[idx_node].idx_node_right = idx_node_right;
        construct_kdtree(
            nodes,
            idx_node_right,
            points,
            idx_point_mid + 1,
            idx_point_end,
            i_depth + 1,
        );
    }
}

/// rasterization of the Kd-tree
fn draw_kdtree(
    pix2rgb: &mut [u8],
    img_width: usize,
    transform_world2pix: &[f32; 9],
    nodes: &[Node],
    idx_node: usize,
    aabb: &[f32; 4],
    i_depth: usize,
) {
    if idx_node >= nodes.len() {
        return;
    }
    let pos = nodes[idx_node].pos;
    if i_depth % 2 == 0 {
        del_canvas::rasterize::line2::draw_dda(
            pix2rgb,
            img_width,
            &[pos.x, aabb[1]],
            &[pos.x, aabb[3]],
            transform_world2pix,
            1,
        );
        draw_kdtree(
            pix2rgb,
            img_width,
            transform_world2pix,
            nodes,
            nodes[idx_node].idx_node_left,
            &[aabb[0], aabb[1], pos.x, aabb[3]],
            i_depth + 1,
        );
        draw_kdtree(
            pix2rgb,
            img_width,
            transform_world2pix,
            nodes,
            nodes[idx_node].idx_node_right,
            &[pos.x, aabb[1], aabb[2], aabb[3]],
            i_depth + 1,
        );
    } else {
        del_canvas::rasterize::line2::draw_dda(
            pix2rgb,
            img_width,
            &[aabb[0], pos.y],
            &[aabb[2], pos.y],
            transform_world2pix,
            1,
        );
        draw_kdtree(
            pix2rgb,
            img_width,
            transform_world2pix,
            nodes,
            nodes[idx_node].idx_node_left,
            &[aabb[0], aabb[1], aabb[2], pos.y],
            i_depth + 1,
        );
        draw_kdtree(
            pix2rgb,
            img_width,
            transform_world2pix,
            nodes,
            nodes[idx_node].idx_node_right,
            &[aabb[0], pos.y, aabb[2], aabb[3]],
            i_depth + 1,
        );
    }
}

/// signed distance from axis-aligned bounding box
/// # Arg
/// * `pos_in` - where the signed distance is evaluated
/// * `aabb` - `[x_min, y_min, x_max, y_max]`
/// # Return
/// signed distance from the `aabb` at `pos_in`  (inside is negative)
fn signed_distance_aabb(pos_in: &Vec2f, aabb: &[f32; 4]) -> f32 {
    let x_center = (aabb[0] + aabb[2]) * 0.5;
    let y_center = (aabb[1] + aabb[3]) * 0.5;
    let x_dist = (pos_in.x - x_center).abs() - (aabb[2] - aabb[0]) * 0.5;
    let y_dist = (pos_in.y - y_center).abs() - (aabb[3] - aabb[1]) * 0.5;
    x_dist.max(y_dist)
}

/// compute nearest position efficiently using kd-tree
/// naive tree traversal approach
/// # Arg
/// * `idx_node_nearest` - the current best nearest node's index
/// * `pos_in` - compute distance from this point
/// * `nodes` - array of nodes
/// * `aabb_node` - Axis-aligned bounding box of this node `[x_min, y_min, x_max, y_max]`
/// * `i_depth` - depth of tree
fn nearest_kdtree_naive(
    idx_node_nearest: &mut usize,
    pos_in: &Vec2f,
    nodes: &[Node],
    idx_node: usize,      // the current node
    aabb_node: &[f32; 4], // aabb of the current node
    i_depth: usize,
) {
    if idx_node >= nodes.len() {
        // this node does not exist
        return;
    }

    // --------------------------------------------------------------
    // write some coe below to cull the branch of the Kd-tree.
    // Check if the region covered this branch does not contain a point that is nearer to the current nearest.
    let dist_to_aabb = signed_distance_aabb(pos_in, aabb_node);
    if dist_to_aabb > 0.0 && dist_to_aabb > (nodes[*idx_node_nearest].pos - pos_in).norm() {
        return;
    }
    // no further edit from here
    // ---------------------------------------------------------------

    let pos_node = nodes[idx_node].pos;

    if (pos_node - pos_in).norm() < (nodes[*idx_node_nearest].pos - pos_in).norm() {
        *idx_node_nearest = idx_node;
    }

    if i_depth % 2 == 0 {
        // division in x direction
        let aabb_west = [aabb_node[0], aabb_node[1], pos_node.x, aabb_node[3]];
        let aabb_east = [pos_node.x, aabb_node[1], aabb_node[2], aabb_node[3]];
        nearest_kdtree_naive(
            idx_node_nearest,
            pos_in,
            nodes,
            nodes[idx_node].idx_node_left,
            &aabb_west,
            i_depth + 1,
        );
        nearest_kdtree_naive(
            idx_node_nearest,
            pos_in,
            nodes,
            nodes[idx_node].idx_node_right,
            &aabb_east,
            i_depth + 1,
        );
    } else {
        // division in y-direction
        let aabb_north = [aabb_node[0], pos_node.y, aabb_node[2], aabb_node[3]];
        let aabb_south = [aabb_node[0], aabb_node[1], aabb_node[2], pos_node.y];
        // division in y-direction
        nearest_kdtree_naive(
            idx_node_nearest,
            pos_in,
            nodes,
            nodes[idx_node].idx_node_left,
            &aabb_south,
            i_depth + 1,
        );
        nearest_kdtree_naive(
            idx_node_nearest,
            pos_in,
            nodes,
            nodes[idx_node].idx_node_right,
            &aabb_north,
            i_depth + 1,
        );
    }
}

/// compute nearest position efficiently using kd-tree
/// # Arg
/// * `idx_node_nearest` - the current best nearest node's index
/// * `pos_in` - compute distance from this point
/// * `nodes` - array of nodes
/// * `aabb_node` - Axis-aligned bounding box of this node `[x_min, y_min, x_max, y_max]`
/// * `i_depth` - depth of tree
fn nearest_kdtree_faster(
    idx_node_nearest: &mut usize,
    pos_in: &Vec2f,
    nodes: &[Node],
    idx_node: usize,      // the current node
    aabb_node: &[f32; 4], // aabb of the current node
    i_depth: usize,
) {
    if idx_node >= nodes.len() {
        // this node does not exist
        return;
    }

    // --------------------------------------------------------------
    // write some coe below to cull the branch of the Kd-tree.
    // Check if the region covered this branch does not contain a point that is nearer to the current nearest.
    let dist_to_aabb = signed_distance_aabb(pos_in, aabb_node);
    if dist_to_aabb > 0.0 && dist_to_aabb > (nodes[*idx_node_nearest].pos - pos_in).norm() {
        return;
    }
    // Write the culling code that is the same as `nearest kdtree naive`

    // no further edit from here
    // ---------------------------------------------------------------

    let pos_node = nodes[idx_node].pos;

    if (pos_node - pos_in).norm() < (nodes[*idx_node_nearest].pos - pos_in).norm() {
        *idx_node_nearest = idx_node;
    }

    // ------------------------------------------------------------------------------------------
    // modify the code below to change the order of the branch evaluation for further acceleration

    if i_depth % 2 == 0 {
        // division in x direction
        let aabb_west = [aabb_node[0], aabb_node[1], pos_node.x, aabb_node[3]];
        let aabb_east = [pos_node.x, aabb_node[1], aabb_node[2], aabb_node[3]];

        if pos_in.x < pos_node.x {
            // Query point is in the left/west branch
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_left,
                &aabb_west,
                i_depth + 1,
            );
            // Only search the right branch if necessary
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_right,
                &aabb_east,
                i_depth + 1,
            );
        } else {
            // Query point is in the right/east branch
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_right,
                &aabb_east,
                i_depth + 1,
            );
            // Only search the left branch if necessary
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_left,
                &aabb_west,
                i_depth + 1,
            );
        }
    } else {
        // division in y-direction
        let aabb_north = [aabb_node[0], pos_node.y, aabb_node[2], aabb_node[3]];
        let aabb_south = [aabb_node[0], aabb_node[1], aabb_node[2], pos_node.y];

        // Choose which branch to visit first based on the query point's position
        if pos_in.y < pos_node.y {
            // Query point is in the bottom/south branch
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_left,
                &aabb_south,
                i_depth + 1,
            );
            // Only search the top branch if necessary
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_right,
                &aabb_north,
                i_depth + 1,
            );
        } else {
            // Query point is in the top/north branch
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_right,
                &aabb_north,
                i_depth + 1,
            );
            // Only search the bottom branch if necessary
            nearest_kdtree_faster(
                idx_node_nearest,
                pos_in,
                nodes,
                nodes[idx_node].idx_node_left,
                &aabb_south,
                i_depth + 1,
            );
        }
    }
    // do not modify from there
    // --------------------------
}

fn make_problem(aabb2: &[f32; 4], num_point: usize) -> Vec<Node> {
    use rand::Rng;
    let mut rng = rand::rng(); // random number generator
    let mut xys: Vec<Vec2f> = vec![Vec2f::new(0., 0.); num_point]; // set number of particles
    for p in &mut xys {
        // // set coordinates
        p.x = aabb2[0] + (aabb2[2] - aabb2[0]) * rng.random::<f32>();
        p.y = aabb2[1] + (aabb2[3] - aabb2[1]) * rng.random::<f32>();
    }
    let mut nodes = Vec::<Node>::with_capacity(xys.len());
    nodes.push(Node::default());
    let nxy = xys.len();
    construct_kdtree(&mut nodes, 0, &mut xys, 0, nxy, 0);
    nodes
}

fn visualize(nodes: &[Node], aabb2: &[f32; 4]) -> anyhow::Result<()> {
    let img_size = 300;
    // mapping from normalized device coordinate (i.e., [-1,+1]^2) to pixel coordinate (i.e., (0,0)x(W,H))
    // This affine transformation matrix is stored in the column major order
    let transform_ndc2pix = [
        img_size as f32 * 0.5,
        0.,
        0.,
        0.,
        -(img_size as f32) * 0.5,
        0.,
        0.5 * (img_size as f32),
        img_size as f32 * 0.5,
        1.,
    ];

    let mut canvas = del_canvas::canvas_gif::Canvas::new(
        "problem0.gif",
        (img_size, img_size),
        &[0xffffff, 0x000000, 0xff0000, 0xff00ff, 0x0000ff],
    )?;
    for i_frame in 0..300 {
        let p1 = {
            let time = i_frame as f32 * 0.01;
            let x = (time * 11.0).sin();
            let y = (time * 7.0).cos();
            Vec2f::new(x, y)
        };
        let mut idx_node_nearest = 0usize;
        nearest_kdtree_naive(&mut idx_node_nearest, &p1, nodes, 0, aabb2, 0);
        canvas.clear(0);
        draw_kdtree(
            &mut canvas.data,
            canvas.width,
            transform_ndc2pix.as_slice().try_into()?,
            nodes,
            0,
            aabb2,
            0,
        );
        // draw bounding box
        del_canvas::rasterize::aabb2::stroke_dda(
            &mut canvas.data,
            canvas.width,
            aabb2,
            transform_ndc2pix.as_slice().try_into()?,
            1,
        );
        del_canvas::rasterize::line2::draw_dda(
            &mut canvas.data,
            canvas.width,
            &[p1.x, p1.y],
            &[nodes[idx_node_nearest].pos.x, nodes[idx_node_nearest].pos.y],
            &transform_ndc2pix,
            3,
        );
        for node in nodes.iter() {
            del_canvas::rasterize::xy::paint_nxn_pixels(
                &mut canvas.data,
                canvas.width,
                &[node.pos[0], node.pos[1]],
                &transform_ndc2pix,
                2,
                3,
            );
        }
        del_canvas::rasterize::xy::paint_nxn_pixels(
            &mut canvas.data,
            canvas.width,
            &[p1.x, p1.y],
            &transform_ndc2pix,
            4,
            3,
        );
        canvas.write();
    }
    Ok(())
}

#[derive(clap::Parser, Debug)]
struct Args {
    #[arg(long, default_value_t = 1000)]
    num_particle: usize,

    #[arg(long, default_value_t = false)]
    vis: bool,
}
fn main() -> anyhow::Result<()> {
    use clap::Parser;
    let args = Args::parse();
    // shape of an axis aligned bounding box (AABB) [x_min, y_min, x_max, y_max].
    let aabb2 = [-0.75f32, -0.75, 0.75, 0.75];
    println!("task04: Accelerated Nearest Search using Kd-Tree");
    println!("number of particles: {}", args.num_particle);
    let nodes = make_problem(&aabb2, args.num_particle);
    if args.vis {
        visualize(&nodes, &aabb2)?;
        return Ok(());
    }
    // -----------------
    use rand::Rng;
    let mut rng = rand::rng(); // random number generator
    let input_points: Vec<Vec2f> = (0..args.num_particle)
        .map(|_| {
            let x = aabb2[0] + (aabb2[2] - aabb2[0]) * rng.random::<f32>();
            let y = aabb2[1] + (aabb2[3] - aabb2[1]) * rng.random::<f32>();
            Vec2f::new(x, y)
        })
        .collect();
    // -------------
    let time_measurement = std::time::Instant::now();
    let mut result_kdtree_naive = vec![];
    for pos_in in &input_points {
        let mut idx_node_nearest = 0usize;
        nearest_kdtree_naive(&mut idx_node_nearest, pos_in, &nodes, 0, &aabb2, 0);
        result_kdtree_naive.push(idx_node_nearest);
    }
    println!(
        "time for kd-tree naive: {} ms",
        time_measurement.elapsed().as_millis()
    );
    // -------------
    let time_measurement = std::time::Instant::now();
    let mut result_kdtree_faster = vec![];
    for pos_in in &input_points {
        let mut idx_node_nearest = 0usize;
        nearest_kdtree_faster(&mut idx_node_nearest, pos_in, &nodes, 0, &aabb2, 0);
        result_kdtree_faster.push(idx_node_nearest);
    }
    println!(
        "time for kd-tree faster: {} ms",
        time_measurement.elapsed().as_millis()
    );
    let time_measurement = std::time::Instant::now();
    // ------------
    let mut result_bruteforce = vec![];
    for pos_in in &input_points {
        let mut idx_node_nearest = 0;
        let mut dist_nearest = (nodes[0].pos - pos_in).norm();
        for (i_node, node) in nodes.iter().enumerate().skip(1) {
            let dist = (node.pos - pos_in).norm();
            if dist < dist_nearest {
                idx_node_nearest = i_node;
                dist_nearest = dist;
            }
        }
        result_bruteforce.push(idx_node_nearest);
    }
    println!(
        "time for brute force: {} ms",
        time_measurement.elapsed().as_millis()
    );
    result_bruteforce
        .iter()
        .zip(result_kdtree_naive.iter())
        .for_each(|(&idx0, &idx1)| {
            assert_eq!(idx0, idx1, "there is bug in your code: {} {}", idx0, idx1)
        });
    result_bruteforce
        .iter()
        .zip(result_kdtree_faster.iter())
        .for_each(|(&idx0, &idx1)| {
            assert_eq!(idx0, idx1, "there is bug in your code: {} {}", idx0, idx1)
        });
    Ok(())
}
