/*!
 * Task01:  Implicit Time Integration
 */

/// Explicit time integration
///
/// # Argument
/// * `rv0` - pair of altitude from the center of the planet and its velocity
/// * `dt` - time step
fn time_integration_explicit(rv0: &mut (f32, f32), dt: f32) {
    let r0 = rv0.0; // current altitude
    let v0 = rv0.1; // current altitude velocity
    let f0 = -1f32 / (r0 * r0); // force based on the **current** altitude
    *rv0 = (r0 + dt * v0, v0 + dt * f0);
}

/// inverse of 2x2 matrix
/// # Argument
/// * `a` : 2x2 matrix with `f32` elements
///
/// # Return
/// Option of unversed 2x2 matrix. The option is `None` is the input is singular.
/// Read https://doc.rust-lang.org/rust-by-example/std/option.html for the details of the Option syntax.
fn inverse_matrix_2x2(a: &[[f32; 2]; 2]) -> Option<[[f32; 2]; 2]> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det == 0.0 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [a[1][1] * inv_det, -a[0][1] * inv_det],
        [-a[1][0] * inv_det, a[0][0] * inv_det],
    ])
}

/// matrix-vector multiplication for 2x2 matrix.
fn mult_mat2_vec(a: &[[f32; 2]; 2], b: &[f32; 2]) -> [f32; 2] {
    [
        a[0][0] * b[0] + a[0][1] * b[1],
        a[1][0] * b[0] + a[1][1] * b[1],
    ]
}

/// Implicit time integration
///
/// # Argument
/// * `rv0` - pair of altitude from the center of the planet and its velocity
/// * `dt` - time step
fn time_integration_implicit(rv0: &mut (f32, f32), dt: f32) {
    let r0 = rv0.0;
    let v0 = rv0.1;
    // ----------------------
    // write some code below

    let dfdr = 2f32 / (r0 * r0 * r0); // hint!

    let a_mat = [[1.0, -dt], [-dfdr * dt, 1.0]];
    let b_vec = [dt * v0, -dt / (r0 * r0)];
    let a_mat_inv = inverse_matrix_2x2(&a_mat).unwrap();
    let res = mult_mat2_vec(&a_mat_inv, &b_vec);
    *rv0 = (r0 + res[0], v0 + res[1]);

    // no further edit from here
    // ----------------------
}

/// When the center of the ball is below the surface of the planet,
/// we set the ball at the surface while setting the velocity such that there is no energy loss
fn reflection(rv0: &mut (f32, f32)) {
    if rv0.0 > 0.5f32 {
        return;
    }
    let r0 = rv0.0;
    let v0 = rv0.1;
    let energy0 = 0.5f32 * v0 * v0 - 1f32 / r0; // energy before reflection
    let r1 = 0.5f32;
    let v1 = (2f32 * energy0 + 4f32).max(0.).sqrt();
    let _energy1 = 0.5f32 * v1 * v1 - 1f32 / r1; // energy before reflection

    // dbg!(_energy1, energy0);

    *rv0 = (r1, v1); // energy preserving reflection
}

fn draw_falling_object(
    pix2coloridx: &mut [u8],
    img_width: usize,
    pos_explicit: &[f32; 2],
    trajectory: &[[f32; 2]],
    transform_ndc2pix: &[f32; 9],
    coloridx: u8,
) {
    del_canvas::rasterize::circle2::stroke_dda(
        pix2coloridx,
        img_width,
        pos_explicit,
        0.05,
        transform_ndc2pix,
        coloridx,
    );
    for i_frame in 0..trajectory.len() - 1 {
        del_canvas::rasterize::line2::draw_dda(
            pix2coloridx,
            img_width,
            &trajectory[i_frame],
            &trajectory[i_frame + 1],
            transform_ndc2pix,
            coloridx,
        );
    }
}

fn main() -> anyhow::Result<()> {
    println!("Task01: Implicit Time Integration");
    let img_size = 300;
    // define canvas for gif animation
    let mut canvas = del_canvas::canvas_gif::Canvas::new(
        "result.gif",
        (img_size, img_size),
        &[
            0xffffff, // color_0
            0x000000, // color_1
            0xff0000, // color_2
            0x0000ff, // color_3
        ],
    )?;
    // mapping from normalized device coordinate (NDC), i.e., [-1,+1]^2 to pixel coordinate, i.e., (0,0)x(W,H).
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
    let mut rv_explicit: (f32, f32) = (0.7, 0.0); // altitude & velocity for explicit time integration
    let mut trajectory_explicit: Vec<[f32; 2]> = vec![];
    let mut rv_implicit: (f32, f32) = (0.7, 0.2); // altitude & velocity for implicit time integration
    let mut trajectory_implicit: Vec<[f32; 2]> = vec![];
    let dt = 0.02;
    let mut time = 0f32;
    for i_frame in 0..1000 {
        //
        time_integration_explicit(&mut rv_explicit, dt);
        reflection(&mut rv_explicit);
        let pos_explicit = [rv_explicit.0 * time.cos(), rv_explicit.0 * time.sin()];
        trajectory_explicit.push(pos_explicit);
        //
        time_integration_implicit(&mut rv_implicit, dt);
        reflection(&mut rv_implicit);
        let pos_implicit = [rv_implicit.0 * time.cos(), rv_implicit.0 * time.sin()];
        trajectory_implicit.push(pos_implicit);
        time += dt;
        // ------------------------------
        // visualization start from here
        if i_frame % 10 == 0 {
            canvas.clear(0); // clear with color_0 (white)
                             // draw planet at the center
            del_canvas::rasterize::circle2::stroke_dda(
                &mut canvas.data,
                canvas.width,
                &[0.0, 0.0],
                0.5,
                &transform_ndc2pix,
                1, // color_1 (black)
            );
            draw_falling_object(
                &mut canvas.data,
                canvas.width,
                &pos_explicit,
                &trajectory_explicit,
                &transform_ndc2pix,
                2,
            ); // color_2 ()
            draw_falling_object(
                &mut canvas.data,
                canvas.width,
                &pos_implicit,
                &trajectory_implicit,
                &transform_ndc2pix,
                3,
            );
            canvas.write();
        }
    }
    Ok(())
}
