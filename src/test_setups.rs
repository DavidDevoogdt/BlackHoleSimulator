use crate::curved_space;
use curved_space::Metric;
use curved_space::SpaceObject;

use crate::python_interface;
use crate::ray_tracer;
use crate::black_body;

extern crate image;
use self::image::{ImageBuffer, RgbImage};
use std::path::Path;


//
pub fn launch_parallel_photons_kerr(){
    let M= 1.0;
    let J = 1.0;
    let a = J;
    let r_s: f64 =2.0;

    let direction = [0.0, 1.0, 0.0];
    let center =[0.0,  0.0, 1.0];
    let inc_vector    = [0.0, 1.0, 0.0]; 


    let center =[-5.0,  1.0, 0.0];
    let inc_vector = [0.0, 0.0, 0.1];


    let accuracy = 1e-5;

    let metric = curved_space::new_kerr_metric(J,  accuracy);
   
    let black_sphere = ray_tracer::Sphere{
        color1: image::Rgb([0,0,0]),
        color2: image::Rgb([0,0,0]),
        radius: 1.0001*(r_s + (r_s.powi(2)-4.0*a.powi(2)).sqrt()  )/2.0,
        divisions: 10.0,
    };

    //let image = image::open("src_files/equirect.png").unwrap().into_rgb();
    let image = image::open("src_files/ESO_-_Milky_Way.jpg").unwrap().into_rgb();
    //let image = image::open("src_files/download.jpeg").unwrap().into_rgb();
    
    let (xres,yres) = image.dimensions();
 
    let skybox = ray_tracer::SkyboxKerr{
        image: image,
        radius: 20.0,
        x_res: xres as i32,
        y_res : yres as i32,
        phi_offset: -90.0/180.0* std::f64::consts::PI,
        a: a,
    };


    let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![
        Box::new( black_sphere),
        Box::new( skybox),

    ];
   
    let mut ray_tracer = ray_tracer::new_parallel( 
        &metric,
        &col_objects,
        direction,
        center ,
        inc_vector,
        100,
        1000,
    );

    ray_tracer.run_simulation( accuracy );
    ray_tracer.plot_paths("xyz");

}

//////
/// 


pub fn ray_trace_kerr(){

    let M= 1.0;
    let J = 0.8;
    let a = J;
    let r_s: f64 =2.0;

    let accuracy = 1e-6;


    let metric = curved_space::new_kerr_metric(J,  accuracy);


    let camera = ray_tracer::Camera{ 
        pos : [-30.0,0.0,5.0],
        direction : [1.0,0.0,-5.0/30.0],
        x_res : 1920/2,
        y_res : 1080/2,
        distance: 0.07,
        width: 0.16,
        height : 0.09,
        rotation_angle :0.0/360.0*(2.0*std::f64::consts::PI),
    };

    //let image = image::open("src_files/milky_way_equirectangular.png").unwrap().into_rgb();
    //let image = image::open("src_files/equirect.png").unwrap().into_rgb();
    let image = image::open("src_files/ESO_-_Milky_Way.jpg").unwrap().into_rgb();
    //let image = image::open("src_files/download.jpeg").unwrap().into_rgb();
    
    let (xres,yres) = image.dimensions();

    //todo: implement proper skybox for kerr metric

    let skybox = ray_tracer::SkyboxKerr{
        image: image,
        radius: 80.0,
        x_res: xres as i32,
        y_res : yres as i32,
        phi_offset: -0.0/180.0* std::f64::consts::PI,
        a: a,
    };
   
    let black_sphere = ray_tracer::Sphere{
        color1: image::Rgb([0,0,0]),
        color2: image::Rgb([0,0,0]),
        radius: 1.0001*(r_s + (r_s.powi(2)-4.0*a.powi(2)).sqrt()  )/2.0,
        divisions: 10.0,
    };

    let accretion_disk = ray_tracer::Annulus{
        color1: image::Rgb([255,0,0]),
        color2: image::Rgb([0,0,255]),
        radius1: 3.0*(r_s + (r_s.powi(2)-4.0*a.powi(2)).sqrt()  )/2.0,
        radius2: 7.0*(r_s + (r_s.powi(2)-4.0*a.powi(2)).sqrt()  )/2.0,
        divisions_angular: 40,
        divisions_radial: 8,
    };

    let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![
        Box::new( black_sphere),
        Box::new(skybox),
        Box::new(accretion_disk),
    ];
   
    let mut ray_tracer = ray_tracer::new( 
        camera,
        &col_objects,
        &metric,
        10000 ,
        false,
    );

    ray_tracer.run_simulation( accuracy );

    ray_tracer.generate_image("src_files/kerr.bmp");
}




//////////////////
/// 

pub fn ray_trace_schwarzshild(){

    let precision = 0.01;
    let steps: i32 = (10000.0/precision) as i32;

    let r_s = 1.0;
    //let metric = curved_space::new_schwarzschild_metric(1.0);
    let metric = curved_space::SchwarzschildMetric{
        r_s:r_s,
        delta: precision,
        max_step: 2.0,
    };

   
    let camera = ray_tracer::Camera{ 
        pos : [0.0,-30.0,5.0],
        direction : [0.0,1.0,-5.0/30.0],
        x_res : 1920,
        y_res : 1080,
        distance: 0.07,
        width: 0.16,
        height : 0.09,
        rotation_angle :0.0/360.0*(2.0*std::f64::consts::PI),
    };

    //let image = image::open("src_files/milky_way_equirectangular.png").unwrap().into_rgb();
    //let image = image::open("src_files/equirect.png").unwrap().into_rgb();
    let image = image::open("src_files/ESO_-_Milky_Way.jpg").unwrap().into_rgb();
    //let image = image::open("src_files/download.jpeg").unwrap().into_rgb();
    
    let (xres,yres) = image.dimensions();

    //todo: implement proper skybox for kerr metric

    let skybox = ray_tracer::Skybox{
        image: image,
        radius: 80.0,
        x_res: xres as i32,
        y_res : yres as i32,
        phi_offset: -0.0/180.0* std::f64::consts::PI,
    };
   
    let black_sphere = ray_tracer::Sphere{
        color1: image::Rgb([0,0,0]),
        color2: image::Rgb([0,0,0]),
        radius: 1.0001*r_s,
        divisions: 10.0,
    };

    let accretion_disk = ray_tracer::Annulus{
        color1: image::Rgb([255,0,0]),
        color2: image::Rgb([0,0,255]),
        radius1: 3.0*r_s,
        radius2: 7.0*r_s,
        divisions_angular: 10,
        divisions_radial: 4,
    };

    let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![
        Box::new( black_sphere),
        Box::new(skybox),
        Box::new(accretion_disk),
    ];

    // let colored_sphere = ray_tracer::Sphere{
    //     color1: image::Rgb([0,255,0]),
    //     color2: image::Rgb([255,255,255]),
    //     radius: 1.001*r_s,
    //     divisions: 5.0,
    // };
    
    // let accretion_disk = ray_tracer::Annulus{
    //     color1: image::Rgb([255,0,0]),
    //     color2: image::Rgb([0,0,255]),
    //     radius1: 3.0*r_s,
    //     radius2: 7.0*r_s,
    //     divisions_angular: 10,
    //     divisions_radial: 5,
    // };
    
    
    // let image = image::open("src_files/ESO_-_Milky_Way.jpg").unwrap().into_rgb();
    // //let image = image::open("src/download.jpeg").unwrap().into_rgb();
    // let (xres,yres) = image.dimensions();

    // let skybox = ray_tracer::Skybox{
    //     image: image,
    //     radius: 20.0,
    //     x_res: xres as i32,
    //     y_res : yres as i32,
    //     phi_offset: std::f64::consts::PI,
    // };
   
    // let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![
    //     //Box::new(black_sphere),
    //     Box::new(colored_sphere),
    //     Box::new(accretion_disk),
    //     Box::new(skybox)
    // ];
   
    let mut ray_tracer = ray_tracer::new( 
        camera,
        &col_objects,
        &metric,
        steps ,
        false,
    );

    
    ray_tracer.run_simulation( 1e1 );

    //ray_tracer.plot_paths();
    ray_tracer.generate_image("src_files/schw800x800.bmp");

    //ray_tracer.plot_paths("xyz");

}


pub fn ray_trace_minkowski(){
    let r_s = 0.1;

    let metric = curved_space::new_minkowski_metric(
        0.01,
    );

    let camera = ray_tracer::Camera{ 
        pos : [ -4.0,-4.0,0.0],
        direction : [1.0,1.0,0.0],
        x_res : 200,
        y_res : 200,
        distance: 0.3,
        height : 0.3,
        width: 0.3,
        rotation_angle :0.0,
     };

    let _colored_sphere = ray_tracer::SphereCart{
        color1: image::Rgb([255,0,0]),
        color2: image::Rgb([0,0,255]),
        radius: 1.05*r_s,
        divisions: 10.0,
    };
    
    let _accretion_disk = ray_tracer::AnnulusCart{
        color1: image::Rgb([255,0,0]),
        color2: image::Rgb([0,0,255]),
        radius1: 3.0*r_s,
        radius2: 6.0*r_s,
        divisions_angular: 10,
        divisions_radial: 5,
    };  

    let image = image::open("src_files/ESO_-_Milky_Way.jpg").unwrap().into_rgb();
    //let image = image::open("src/download.jpeg").unwrap().into_rgb();
    let (xres,yres) = image.dimensions();

    let skybox = ray_tracer::SkyboxCart{
        image: image,
        radius: 10.0,
        x_res: xres as i32,
        y_res : yres as i32,
        phi_offset: 90.0/360.0 *std::f64::consts::PI,
    };
   
        
    //let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![Box::new(colored_sphere), Box::new(accretion_disk), Box::new(skybox) ];
    //let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![Box::new(colored_sphere), Box::new(accretion_disk) ];
    let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![Box::new(skybox) ];

    let mut ray_tracer = ray_tracer::new( 
        camera,
        &col_objects,
        &metric,
        2000,
        false
    );

    
    ray_tracer.run_simulation( 1e-3  );


    ray_tracer.generate_image("files/flat800x800.bmp");

    //ray_tracer.plot_paths();

}


// blackbodytester

pub fn generate_blackbody(temp_min: f64, temp_max:f64, steps: u32) {
    let cie = black_body::CieLookup::new();
    

    let ysteps = ( (steps as f64) /10.0).ceil() as u32;

    let mut img: RgbImage = ImageBuffer::new(steps, ysteps);

    for i in 0..steps {
        let t = temp_min + (temp_max-temp_min)*(i as f64)/(steps as f64);

        let col = cie.get_black_body_intensity(t, -200.0);
        for j in 0..ysteps {
            img.put_pixel(i, j , col );
        } 
        
    }

    let p = Path::new("files/blackbody.bmp");

    img.save_with_format( p,  image::ImageFormat::Bmp).unwrap();

}

#[test]
fn test_wavelentgh_convo (){
    let cie = black_body::CieLookup::new();

    let arr : [f64;3] = [435.8,546.1,700.0];

    for (_,x) in arr.iter().enumerate() {
        let col = cie.wavelength_to_rgb(*x);
        println!("lambda: {} rgb: {},{},{}",x ,col[0], col[1],col[2] );
    }

 }


 pub fn photon_orbit(){

    let r_s = 0.1;
    //let metric = curved_space::new_schwarzschild_metric(1.0);
    let metric = curved_space::SchwarzschildMetric{
        r_s:r_s,
        delta: 1.0/32.0,
        max_step: 2.0,
    };


    let camera = ray_tracer::Camera{ 
        pos : [ -1.5*r_s,0.0,0.0],
        direction : [0.0,0.0,1.0],
        x_res : 1,
        y_res : 1,
        distance: 0.05,
        width: 0.01,
        height : 0.1,
        rotation_angle :0.0/360.0*(2.0*std::f64::consts::PI),
    };

    let col_objects : Vec< Box< dyn ray_tracer::CollsionObject> > = vec![];

    let mut ray_tracer = ray_tracer::new( 
        camera,
        &col_objects,
        &metric,
        20000 ,
        true,
    );

    ray_tracer.run_simulation( 1e-3 );


    ray_tracer.plot_paths( "xz");

 }


pub fn test_coordinate_system_2(){

    let metric = curved_space::new_kerr_metric(1.0, 1e-1);
 
    let  [coor,mom] = [
        [0.0,-1.3,2.5,0.4],
        [1.0,0.0,1.0,0.0]
    ];

    
    let phot = metric.spawn_space_object_from_cartesian( coor, mom,0.0 );
    let [coor2,mom2] = phot.calculate_cartesian_coordinates_and_momenta();

    println!("[{:.4},{:.4},{:.4},{:.4}]->[{:.4},{:.4},{:.4},{:.4}]", coor[0],coor[1],coor[2],coor[3], coor2[0],coor2[1],coor2[2],coor2[3] );
    println!("[{:.4},{:.4},{:.4},{:.4}]->[{:.4},{:.4},{:.4},{:.4}]", mom[0],mom[1],mom[2],mom[3], mom2[0],mom2[1],mom2[2],mom2[3] );
       

}

pub fn test_coordinate_system_3(){

    let  [coor,mom] = [
        [-1.3,2.31,0.4],
        [1.0,3.0,-0.68]
    ];

    let [coor2,mom2] =curved_space::cart_to_kerr( &coor,&mom,0.5 );
    let [coor3,mom3] =curved_space::kerr_to_cart( &coor2,&mom2, 0.5 );

    println!("[{:.4},{:.4},{:.4}]->[{:.4},{:.4},{:.4}]", coor[0],coor[1],coor[2], coor3[0],coor3[1],coor3[2] );
    println!("[{:.4},{:.4},{:.4}]->[{:.4},{:.4},{:.4}]", mom[0],mom[1],mom[2], mom3[0],mom3[1],mom3[2] );
       

}

pub fn test_coordinate_system(){

    let metr = curved_space::new_kerr_metric(0.32, 0.0);
 
    let  [coor,mom] = [
        [1.0,-1.28,2.85,0.4],
        [1.0,1.0,3.1,-0.68]
    ];

    let co_mom = metr.to_covariant_vector(&coor, &mom);
    let contra_mom = metr.to_contravariant_vector(&coor, &co_mom);


    println!("[{:.4},{:.4},{:.4}{:.4}]->[{:.4},{:.4},{:.4},{:.4}]", mom[0],mom[1],mom[2],mom[3], contra_mom[0],contra_mom[1],contra_mom[2],contra_mom[3] );
       

}
