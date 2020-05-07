use crate::curved_space;
use crate::python_interface;
use crate::ray_tracer;

extern crate chrono;
use self::chrono::Local;

use std::fs;

extern crate image;

#[allow(non_snake_case)]
pub fn make_kerr_video() {
    // setup properties of metric
    //let M= 1.0;
    let J = 0.7;
    let a = J;
    let r_s: f64 = 2.0;

    let accuracy = 1e-7;
    let metric = curved_space::new_kerr_metric(J, accuracy);

    // unchanged camera properties
    let x_res = 1920;
    let y_res = 1080;
    let distance = 0.07;
    let width = 0.16;
    let height = 0.09;
    let rotation_angle = 0.0 / 360.0 * (2.0 * std::f64::consts::PI);

    // setup environment

    //let image = image::open("src_files/milky_way_equirectangular.png").unwrap().into_rgb();
    //let image = image::open("src_files/equirect.png").unwrap().into_rgb();
    let image = image::open("src_files/ESO_-_Milky_Way.jpg")
        .unwrap()
        .into_rgb();
    //let image = image::open("src_files/download.jpeg").unwrap().into_rgb();

    let (xres, yres) = image.dimensions();

    //todo: implement proper skybox for kerr metric

    let skybox = ray_tracer::SkyboxKerr {
        image: image,
        radius: 80.0,
        x_res: xres as i32,
        y_res: yres as i32,
        phi_offset: -0.0 / 180.0 * std::f64::consts::PI,
        a: a,
    };

    let black_sphere = ray_tracer::Sphere {
        color1: image::Rgb([0, 0, 0]),
        color2: image::Rgb([0, 0, 0]),
        radius: 1.0001 * (r_s + (r_s.powi(2) - 4.0 * a.powi(2)).sqrt()) / 2.0,
        divisions: 10.0,
    };

    let accretion_disk = ray_tracer::Annulus {
        color1: image::Rgb([255, 0, 0]),
        color2: image::Rgb([0, 0, 255]),
        radius1: 3.0 * (r_s + (r_s.powi(2) - 4.0 * a.powi(2)).sqrt()) / 2.0,
        radius2: 7.0 * (r_s + (r_s.powi(2) - 4.0 * a.powi(2)).sqrt()) / 2.0,
        divisions_angular: 40,
        divisions_radial: 8,
    };

    let col_objects: Vec<Box<dyn ray_tracer::CollsionObject>> = vec![
        Box::new(black_sphere),
        Box::new(skybox),
        Box::new(accretion_disk),
    ];

    // trajectory properties:
    let r_trajectory: f64 = 30.0;

    let start = -10.0 / 180.0 * std::f64::consts::PI;
    let stop = 20.0 / 180.0 * std::f64::consts::PI;
    let step = 0.05 / 180.0 * std::f64::consts::PI;

    let num = ((stop - start) / step) as i32;

    let base_name = "kerr3";

    let start_time = Local::now();

    fs::create_dir_all("video_files").unwrap();

    for n in 0..num {
        let phi = start + (n as f64) * step;

        let camera = ray_tracer::Camera {
            x_res,
            y_res,
            distance,
            width,
            height,
            rotation_angle,
            pos: [-r_trajectory * phi.cos(), -r_trajectory * phi.sin(), 3.0],
            direction: [phi.cos(), phi.sin(), -3.0 / r_trajectory],
        };

        let mut ray_tracer = ray_tracer::new(camera, &col_objects, &metric, 20000, false);

        ray_tracer.run_simulation(accuracy);
        ray_tracer.generate_image(&format!("video_files/{}{:04}.bmp", base_name, n));

        let duration = Local::now().signed_duration_since(start_time);
        let total_duration = duration / (n + 1) * num;
        let finish_time = start_time + total_duration;

        let days = total_duration.num_days();
        let hours = total_duration.num_hours() - 24 * days;
        let minutes = total_duration.num_minutes() - 60 * total_duration.num_hours();
        let seconds = total_duration.num_seconds() - 60 * total_duration.num_minutes();

        print!("{:.03}%: {} of {}, duration: {} day(s) {:02}:{:02}:{:02}, estimated finish time: {} \n", ((n+1) as f64)/(num as f64)*100.0 ,n+1,num,days,hours,minutes,seconds,finish_time.format("%Y-%m-%d %H:%M:%S"));
    }

    python_interface::make_video(num, base_name, 24);
}
