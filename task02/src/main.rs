//! task02: Linear Momentum conservation

type Vec2f = nalgebra::Vector2<f32>;

/// particle class (radius is zero)
/// * `pos` - 2D position
/// * `velo` - 2D velocity
#[derive(Clone)]
struct Particle {
    pos: Vec2f,
    velo: Vec2f,
}

impl Particle {
    fn default() -> Self {
        Particle {
            pos: nalgebra::zero(),
            velo: nalgebra::zero(),
        }
    }
}

///  collision between a circle and plane
fn collision_circle_plane(
    pos: &mut Vec2f,
    velo: &mut Vec2f,
    rad: f32,
    plane_org: &Vec2f,
    plane_nrm: &Vec2f,
) {
    let height = (*pos - plane_org).dot(plane_nrm) - rad;
    if height > 0. {
        return;
    }
    *pos -= height * 2. * plane_nrm;
    let velo_perp = velo.dot(plane_nrm);
    *velo -= 2. * velo_perp * plane_nrm;
}

fn collide_particle_ball(
    p: &mut Particle,
    particle_mass: &f32,
    ball_pos: &Vec2f,
    ball_velo: &mut Vec2f,
    ball_mass: f32,
    ball_rad: f32,
) {
    if (p.pos - ball_pos).norm() > ball_rad {
        return;
    }
    let plane_norm = (p.pos - ball_pos).normalize();
    let plane_org = ball_pos + plane_norm * ball_rad;
    let height = (p.pos - plane_org).dot(&plane_norm);
    p.pos -= height * 2. * plane_norm;
    ////////////////////////
    // write some code

    // comment out the line below.
    // p.velo -= 2f32 * (p.velo - *ball_velo).dot(&plane_norm) * plane_norm;

    // References: https://en.wikipedia.org/wiki/Elastic_collision
    *ball_velo += -2f32 * particle_mass / (ball_mass + particle_mass)
        * ball_velo.dot(&plane_norm)
        * plane_norm
        + 2f32 * particle_mass / (ball_mass + particle_mass) * p.velo.dot(&plane_norm) * plane_norm;
    p.velo += -2f32 * ball_mass / (ball_mass + particle_mass)
        * p.velo.dot(&plane_norm)
        * plane_norm
        + 2f32 * ball_mass / (ball_mass + particle_mass) * ball_velo.dot(&plane_norm) * plane_norm;

    // no edit from here
    ////////////////////////
}

///  collision between a circle and box boundary
fn collision_against_aabb2(pos: &mut Vec2f, radius: f32, velo: &mut Vec2f, aabb2: &[f32; 4]) {
    // collision against right wall
    collision_circle_plane(
        pos,
        velo,
        radius,
        &Vec2f::new(aabb2[0], 0.),
        &Vec2f::new(1., 0.),
    );
    // collision against left wall
    collision_circle_plane(
        pos,
        velo,
        radius,
        &Vec2f::new(aabb2[2], 0.),
        &Vec2f::new(-1., 0.),
    );
    // collision against top wall
    collision_circle_plane(
        pos,
        velo,
        radius,
        &Vec2f::new(0., aabb2[1]),
        &Vec2f::new(0., 1.),
    );
    // collision against bottom wall
    collision_circle_plane(
        pos,
        velo,
        radius,
        &Vec2f::new(0., aabb2[3]),
        &Vec2f::new(0., -1.),
    );
}

fn main() -> anyhow::Result<()> {
    println!("task02: Linear Momentum Conservation");
    // shape of a axis aligned bounding box (AABB) [x_min, y_min, x_max, y_max].
    let aabb2 = [-0.75f32, -0.75, 0.75, 0.75];

    let ball_rad: f32 = 0.2;
    let ball_mass: f32 = 10.;
    let mut ball_pos = Vec2f::new(0., 0.);
    let mut ball_velo = Vec2f::new(0., 0.);

    let mut particles: Vec<Particle> = vec![Particle::default(); 100];
    let particle_mass = 1f32;
    use rand::Rng;
    let mut rng = rand::rng(); // random number generator
    for particle in &mut particles {
        particle.pos = [
            aabb2[0] + (aabb2[2] - aabb2[0]) * rng.random::<f32>(),
            aabb2[1] + (aabb2[3] - aabb2[1]) * rng.random::<f32>(),
        ]
        .into();
        particle.velo = [
            rng.random::<f32>() * 2f32 - 1f32,
            rng.random::<f32>() * 2f32 - 1f32,
        ]
        .into();
        particle.velo = particle.velo.normalize();
    }
    let img_size = 300;
    // palette (0:white, 1: black, 2:red, 3: green, 4: blue)
    let mut canvas = del_canvas::canvas_gif::Canvas::new(
        "output.gif",
        (img_size, img_size),
        &[0xffffff, 0x000000, 0xff0000, 0x00ff00, 0x0000ff],
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

    let mut ball_trajectory: Vec<[f32; 2]> = vec![];
    let dt = 0.05;
    for _i_frame in 0..200 {
        ball_pos += ball_velo * dt; // step time for ball
        collision_against_aabb2(&mut ball_pos, ball_rad, &mut ball_velo, &aabb2);
        for particle in &mut particles {
            particle.pos += particle.velo * dt;
            collision_against_aabb2(&mut particle.pos, 0f32, &mut particle.velo, &aabb2);
            collide_particle_ball(
                particle,
                &particle_mass,
                &ball_pos,
                &mut ball_velo,
                ball_mass,
                ball_rad,
            );
        }
        ball_trajectory.push(*ball_pos.as_mut());
        // --------------------
        // below: visualization
        canvas.clear(0);
        // draw boundary box
        del_canvas::rasterize::aabb2::stroke_dda(
            &mut canvas.data,
            canvas.width,
            &aabb2,
            &transform_ndc2pix,
            1, // color_1(black)
        );
        // draw ball
        del_canvas::rasterize::circle2::stroke_dda(
            &mut canvas.data,
            canvas.width,
            &[ball_pos.x, ball_pos.y],
            ball_rad,
            &transform_ndc2pix,
            3, // color_3(green)
        );
        // draw particles
        for particle in &particles {
            del_canvas::rasterize::xy::paint_nxn_pixels(
                &mut canvas.data,
                canvas.width,
                &[particle.pos.x, particle.pos.y],
                &transform_ndc2pix,
                1, // color_1(black)
                2,
            );
        }
        // draw trajectory
        del_canvas::rasterize::polyline2::stroke_dda(
            &mut canvas.data,
            canvas.width,
            &ball_trajectory,
            &transform_ndc2pix,
            4,
        );
        canvas.write();
    }

    Ok(())
}
