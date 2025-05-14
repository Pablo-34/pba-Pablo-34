type Vec2f = nalgebra::Vector2<f32>;

#[derive(Default, Clone)]
struct Particle {
    pos: Vec2f,
    velocity: Vec2f,
    force: Vec2f,
}

/// Gravitational force with softening
/// d relative position
fn gravitational_force(d: &Vec2f) -> Vec2f {
    let eps = 2.0e-3; // softening coefficient
    let r: f32 = (d.norm_squared() + eps * eps).sqrt();
    d.scale(1f32 / (r * r * r))
}

/// For each particle, set summation of gravitational forces
/// from all the other particles in a brute-force way O(N^2)
fn set_force_bruteforce(particles: &mut [Particle]) {
    for ip in 0..particles.len() {
        particles[ip].force = Vec2f::new(0., 0.);
        for jp in 0..particles.len() {
            if ip == jp {
                continue;
            }
            particles[ip].force += gravitational_force(&(particles[jp].pos - particles[ip].pos));
        }
    }
}

/// position to grid coordinate
/// # Args
/// * `pos` - input position
/// * `box_size` - size of square box
/// * `num_div` - number of division for grid
/// # Return
/// cell index
fn cell_index_from_position(pos: &Vec2f, box_size: f32, num_div: usize) -> usize {
    let h_inv = num_div as f32 / box_size;
    let ix = ((pos.x + box_size * 0.5f32) * h_inv).floor() as isize;
    let iy = ((pos.y + box_size * 0.5f32) * h_inv).floor() as isize;
    if ix < 0 || ix >= num_div as isize {
        // outside grid
        return usize::MAX;
    };
    if iy < 0 || iy >= num_div as isize {
        // outside grid
        return usize::MAX;
    };
    (iy * num_div as isize + ix) as usize
}

/// data structure for acceleration
struct Acceleration {
    box_size: f32,  // size of the box
    num_div: usize, // number of division of each side of the box

    // `idx2ipic` and `cell2idx` are jagged array data structure,
    // storing pairs (particle index, cell index) for each cell
    cell2idx: Vec<usize>,          // index of jagged array
    idx2ipic: Vec<(usize, usize)>, // data of jagged array

    cell2cg: Vec<Vec2f>, // the center of gravity of each grid
}

impl Acceleration {
    fn construct(&mut self, particles: &[Particle]) {
        let num_cell = self.num_div * self.num_div;
        self.idx2ipic.reserve(particles.len());
        self.idx2ipic.clear(); // empty vector while retaining memory
        for (i_particle, particle) in particles.iter().enumerate() {
            let i_grid = cell_index_from_position(&particle.pos, self.box_size, self.num_div);
            self.idx2ipic.push((i_particle, i_grid));
        }
        self.idx2ipic.sort_by(|&a, &b| a.1.cmp(&b.1)); // sort by the cell index
                                                       // make jagged array structure
        self.cell2idx.resize(num_cell + 1, 0);
        self.cell2idx.fill(0);
        for &(_i_particle, i_cell) in &self.idx2ipic {
            if i_cell == usize::MAX {
                // this particle is outside grid
                continue;
            }
            // count the number of particle in the cell
            self.cell2idx[i_cell + 1] += 1;
        }
        for i_grid in 0..num_cell {
            // compute prefix sum
            self.cell2idx[i_grid + 1] += self.cell2idx[i_grid];
        }
        // compute the center of the gravity for each grid cell
        self.cell2cg.resize(num_cell, Vec2f::new(0., 0.));
        self.cell2cg.fill(Vec2f::new(0., 0.));
        for i_grid in 0..num_cell {
            for idx in self.cell2idx[i_grid]..self.cell2idx[i_grid + 1] {
                let (i_particle, i_grid) = self.idx2ipic[idx];
                self.cell2cg[i_grid] += particles[i_particle].pos;
            }
            let num_particle_in_cell = self.cell2idx[i_grid + 1] - self.cell2idx[i_grid];
            if num_particle_in_cell == 0 {
                continue;
            }
            self.cell2cg[i_grid] /= num_particle_in_cell as f32;
        }
    }
}

/// For each particle, set summation of gravitational forces
/// from all the other particles in an accelerated way
/// # Args
/// * `particles` - particles
/// * `acc` - acceleration data structure
fn set_force_accelerated(particles: &mut [Particle], acc: &Acceleration) {
    for i_particle in 0..particles.len() {
        // cell index of this particle
        let i_cell =
            cell_index_from_position(&particles[i_particle].pos, acc.box_size, acc.num_div);
        // grid coordinate
        let (ix, iy) = (i_cell % acc.num_div, i_cell / acc.num_div);
        particles[i_particle].force = Vec2f::new(0., 0.);
        // loop over all the grid set force to the particle with index `ip`
        for j_cell in 0..acc.num_div * acc.num_div {
            // grid coordinate of `j_cell`
            let (jx, jy) = (j_cell % acc.num_div, j_cell / acc.num_div);
            if ix.abs_diff(jx) <= 1 && iy.abs_diff(jy) <= 1 {
                // this grid is near to the particle as the grid indexes are close
                for jdx in acc.cell2idx[j_cell]..acc.cell2idx[j_cell + 1] {
                    let j_particle = acc.idx2ipic[jdx].0; // particle index in this cell
                    if i_particle == j_particle {
                        continue;
                    }
                    let diff = particles[j_particle].pos - particles[i_particle].pos;
                    particles[i_particle].force += gravitational_force(&diff);
                }
            } else {
                // far field approximation
                // write the code to approximate the force from particles in this cell.
                let cg = acc.cell2cg[j_cell]; // hint: center of the gravity in this cell
                let num_particle_in_cell = acc.cell2idx[j_cell + 1] - acc.cell2idx[j_cell]; // hint number of the particle in this cell
                particles[i_particle].force +=
                    gravitational_force(&(cg - particles[i_particle].pos))
                        .scale(num_particle_in_cell as f32);
            }
        }
    }
}

/// the arguments of command line
#[derive(clap::Parser, Debug)]
struct Args {
    #[arg(long, default_value_t = 1000)]
    num_particle: usize,

    #[arg(long, default_value_t = false)]
    accelerate: bool,
}

fn main() -> anyhow::Result<()> {
    use clap::Parser;
    let args = Args::parse();
    let num_div = 32;
    println!(
        "Accelerate: {}, Number of Particle: {}",
        args.accelerate, args.num_particle
    );
    let box_size = 1.5f32;
    let mut acc = Acceleration {
        box_size,
        num_div,
        idx2ipic: vec![],
        cell2idx: vec![],
        cell2cg: vec![],
    };
    // shape of a axis aligned bounding box (AABB) [x_min, y_min, x_max, y_max].
    let aabb2 = [
        -box_size * 0.5f32,
        -box_size * 0.5f32,
        box_size * 0.5f32,
        box_size * 0.5f32,
    ];
    let mut particles: Vec<Particle> = vec![Particle::default(); args.num_particle];
    let mut rng = rand::rng(); // random number generator
    for p in &mut particles {
        use rand::Rng;
        p.pos = Vec2f::new(
            aabb2[0] + (aabb2[2] - aabb2[0]) * rng.random::<f32>(),
            aabb2[1] + (aabb2[3] - aabb2[1]) * rng.random::<f32>(),
        );
        let center = [(aabb2[0] + aabb2[2]) * 0.5, (aabb2[1] + aabb2[3]) * 0.5];
        p.velocity = Vec2f::new(100. * (p.pos.y - center[1]), -100. * (p.pos.x - center[0]));
    }
    let img_size = 300;
    let mut canvas = del_canvas::canvas_gif::Canvas::new(
        "output.gif",
        (img_size, img_size),
        &[0xffffff, 0x000000, 0xff0000, 0x00ff00],
    )?;
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

    let start = std::time::Instant::now();
    let dt = 0.00002f32; // time step
    for i_frame in 0..1000 {
        // switch brute-force/accelerated computation here by uncomment/comment below
        if args.accelerate {
            acc.construct(&particles);
            set_force_accelerated(&mut particles, &acc);
        } else {
            set_force_bruteforce(&mut particles);
        }
        for p in &mut particles {
            // leap frog time integration
            p.velocity += p.force * dt; // update velocity
            p.pos += p.velocity * dt; // update position
        }
        if i_frame % 20 != 0 {
            continue;
        }
        // --------------------
        // below: visualization
        canvas.clear(0);
        del_canvas::rasterize::aabb2::stroke_dda(
            &mut canvas.data,
            canvas.width,
            &aabb2,
            &transform_ndc2pix,
            1,
        );
        for particle in &particles {
            del_canvas::rasterize::xy::paint_one_pixel(
                &mut canvas.data,
                canvas.width,
                &[particle.pos.x, particle.pos.y],
                &transform_ndc2pix,
                1,
            );
        }
        canvas.write();
    }
    println!("computation time: {:.2?}", start.elapsed());

    Ok(())
}
