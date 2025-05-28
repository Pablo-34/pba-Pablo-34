type Vec2f = nalgebra::Vector2<f32>;

/// values are defined at the center of pixels.
/// given position `xy`, consider a square connecting the centers of the adjacent pixels.
/// this function returns the index of pixel that is north-west (left up) corner of the square.
/// return None if there is no square connecting center that include position `xy`.
pub fn pixel_north_west_to_xy(xy: &[f32; 2], img_resolution: usize) -> Option<usize> {
    let x = (xy[0] * img_resolution as f32) - 0.5f32;
    let y = (xy[1] * img_resolution as f32) - 0.5f32;

    // ------------------
    // implement some code below.
    if x < 0.0 || y < 0.0 || x >= (img_resolution - 1) as f32 || y >= (img_resolution - 1) as f32 {
        return None;
    }

    let i = x.floor() as usize; // column index
    let j = ((img_resolution - 1) as f32 - y).floor() as usize; // row index

    let index = j * img_resolution + i;

    Some(index)

    // no edit from here
    // -----------------
}

/// test for `pixel_north_west_to_xy` function.
#[test]
fn test_pixel_north_west_to_xy() {
    let ns = [4, 5, 7];
    for &n in &ns {
        let r = 1f32 / n as f32;
        assert_eq!(pixel_north_west_to_xy(&[r, 1.0 - r], n), Some(0));
        assert_eq!(
            pixel_north_west_to_xy(&[r * 0.8, 1.0 - r * 0.8], n),
            Some(0)
        );
        assert_eq!(
            pixel_north_west_to_xy(&[r * 1.2, 1.0 - r * 1.2], n),
            Some(0)
        );
        assert_eq!(pixel_north_west_to_xy(&[1.0 - r, 1.0 - r], n), Some(n - 2));
        assert_eq!(pixel_north_west_to_xy(&[r, 1.0 - 2.0 * r], n), Some(n));
        assert_eq!(pixel_north_west_to_xy(&[r, r], n), Some((n - 2) * n));
        assert_eq!(pixel_north_west_to_xy(&[0.3 * r, r], n), None);
        assert_eq!(pixel_north_west_to_xy(&[r, 0.3 * r], n), None);
        assert_eq!(pixel_north_west_to_xy(&[1.0 - 0.3 * r, r], n), None);
    }
}

/// compute gradient at `xy` when the values are defined at the pixel center.
pub fn gradient(xy: &[f32; 2], img_resolution: usize, pix2val: &[f32]) -> [f32; 2] {
    let Some(i_pix_nw) = pixel_north_west_to_xy(xy, img_resolution) else {
        return [0f32; 2];
    };
    let r = 1.0 / img_resolution as f32;
    let i_pix_sw = i_pix_nw + img_resolution;
    let i_pix_ne = i_pix_nw + 1;
    let i_pix_se = i_pix_nw + img_resolution + 1;
    let pixel_center_sw = [
        ((i_pix_sw % img_resolution) as f32 + 0.5) * r,
        1.0 - ((i_pix_sw / img_resolution) as f32 + 0.5) * r,
    ];
    let rx = (xy[0] - pixel_center_sw[0]) * img_resolution as f32; // horizontal ratio inside square
    let ry = (xy[1] - pixel_center_sw[1]) * img_resolution as f32; // vertical ratio inside square
    let val_nw = pix2val[i_pix_nw];
    let val_ne = pix2val[i_pix_ne];
    let val_sw = pix2val[i_pix_sw];
    let val_se = pix2val[i_pix_se];
    // ---------------------
    // write some code below to compute gradient
    // First, interpolate values along the edges
    let val_north = val_nw * (1.0 - rx) + val_ne * rx; // top edge
    let val_south = val_sw * (1.0 - rx) + val_se * rx; // bottom edge
    let val_west = val_sw * (1.0 - ry) + val_nw * ry; // left edge
    let val_east = val_se * (1.0 - ry) + val_ne * ry; // right edge

    // Compute gradient components
    let grad_x = (val_east - val_west) * img_resolution as f32;
    let grad_y = (val_north - val_south) * img_resolution as f32;

    [grad_x, grad_y]
    // no edit from here
    // -----------------
}

pub fn palette() -> Vec<i32> {
    let mut palette = vec![];
    for i in 0..=254 {
        let v = i as f32 / 254.0;
        let color = del_canvas::colormap::apply_colormap(
            v,
            0.0f32,
            1.0f32,
            del_canvas::colormap::COLORMAP_HOT,
        );
        let c = del_canvas::color::i32_from_f32rgb(color[0], color[1], color[2]);
        palette.push(c);
    }
    {
        // palette color 255 is reserved for point color
        let c = del_canvas::color::i32_from_f32rgb(0., 1., 1.);
        palette.push(c);
    }
    palette
}

fn solve_laplace_gauss_seidel_on_grid(
    pix2val: &mut [f32],
    img_resolution: usize,
    pix2isfix: &[u8],
) {
    assert_eq!(pix2val.len(), img_resolution * img_resolution);
    for i in 0..img_resolution {
        for j in 0..img_resolution {
            let i_pix_center = j * img_resolution + i;
            if pix2isfix[i_pix_center] == 1 {
                continue;
            }
            let val_north = pix2val[j * img_resolution + (i + 1)];
            let val_south = pix2val[j * img_resolution + (i - 1)];
            let val_west = pix2val[(j - 1) * img_resolution + i];
            let val_east = pix2val[(j + 1) * img_resolution + i];
            // ------------------------
            // write some code below

            pix2val[i_pix_center] = (val_north + val_south + val_west + val_east) / 4.0;

            // no edit from here
            // -------------------------------
        }
    }
}

fn dirichlet_energy(pix2val: &[f32], img_resolution: usize) -> f32 {
    let mut eng = 0f32;
    // energy on horizontal edge
    for ix0 in 0..img_resolution - 1 {
        let ix1 = ix0 + 1;
        for iy in 0..img_resolution {
            let val0 = pix2val[iy * img_resolution + ix0];
            let val1 = pix2val[iy * img_resolution + ix1];
            eng += (val0 - val1) * (val0 - val1);
        }
    }
    // energy on vertical edge
    for iy0 in 0..img_resolution - 1 {
        let iy1 = iy0 + 1;
        for ix in 0..img_resolution {
            let val0 = pix2val[iy0 * img_resolution + ix];
            let val1 = pix2val[iy1 * img_resolution + ix];
            eng += (val0 - val1) * (val0 - val1);
        }
    }
    eng
}

#[derive(clap::Parser, Debug)]
struct Args {
    #[arg(long, default_value_t = 256)]
    img_resolution: usize,
}

fn main() -> anyhow::Result<()> {
    use clap::Parser;
    let args = Args::parse();
    println!("task05");
    let img_resolution = args.img_resolution;
    // random initial values at the pixels
    let mut pix2val: Vec<f32> = {
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0);
        (0..img_resolution * img_resolution)
            .map(|_| rng.random::<f32>())
            .collect()
    };
    // set up fixed boundary condition
    let pix2isfix = {
        let mut pix2isfix = vec![0; img_resolution * img_resolution];
        // set fixed boundary condition at the edge of square
        for i in 0..img_resolution {
            let i0 = i * img_resolution + img_resolution - 1;
            let i1 = i * img_resolution;
            let i2 = i;
            let i3 = (img_resolution - 1) * img_resolution + i;
            pix2isfix[i0] = 1;
            pix2isfix[i1] = 1;
            pix2isfix[i2] = 1;
            pix2isfix[i3] = 1;
            pix2val[i0] = 1.;
            pix2val[i1] = 1.;
            pix2val[i2] = 1.;
            pix2val[i3] = 1.;
        }
        // put boundary condition at the center circle
        for ix in 0..img_resolution {
            for iy in 0..img_resolution {
                let x = (ix as f32 + 0.5) / img_resolution as f32;
                let y = 1.0 - (iy as f32 + 0.5) / img_resolution as f32;
                if (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.1 * 0.1 {
                    let ipix = iy * img_resolution + ix;
                    pix2isfix[ipix] = 1;
                    pix2val[ipix] = 0.;
                }
            }
        }
        pix2isfix
    };
    // transformation unit square [0,1]^2 to [0,0,img_size, img_size]
    let transform_usq2pix = [
        img_resolution as f32,
        0.,
        0.,
        0.,
        -(img_resolution as f32),
        0.,
        0.,
        img_resolution as f32,
        1.,
    ];
    let mut canvas = del_canvas::canvas_gif::Canvas::new(
        "output.gif",
        (img_resolution, img_resolution),
        &palette(),
    )?;

    for i_iteration in 0..1000 {
        if i_iteration % 10 == 0 {
            canvas.clear(0);
            for (i_pix, &val) in pix2val.iter().enumerate() {
                canvas.data[i_pix] = (val * 254.0) as u8;
            }
            canvas.write();
        }
        let dirichlet_eng = dirichlet_energy(&pix2val, img_resolution);
        println!("{} {}", i_iteration, dirichlet_eng);
        solve_laplace_gauss_seidel_on_grid(&mut pix2val, img_resolution, &pix2isfix);
    }
    let mut point2xy: Vec<Vec2f> = {
        use rand::Rng;
        use rand::SeedableRng;
        let mut rng = rand_chacha::ChaChaRng::seed_from_u64(0);
        (0..1000)
            .map(|_| {
                let x: f32 = rng.random();
                let y: f32 = rng.random();
                Vec2f::new(x, y)
            })
            .collect()
    };
    for _i_frame in 0..100 {
        for point in &mut point2xy {
            let grad = gradient(&[point.x, point.y], img_resolution, &pix2val);
            // move the particles negative direction of the gradient (heat flow)
            point.x -= grad[0] * 0.001;
            point.y -= grad[1] * 0.001;
        }
        canvas.clear(0);
        for (i_pix, &val) in pix2val.iter().enumerate() {
            canvas.data[i_pix] = (val * 254.0) as u8;
        }
        for xy in &point2xy {
            del_canvas::rasterize::xy::paint_nxn_pixels(
                &mut canvas.data,
                canvas.width,
                &[xy.x, xy.y],
                &transform_usq2pix,
                255, // blue
                2,
            );
        }
        canvas.write();
    }
    Ok(())
}
