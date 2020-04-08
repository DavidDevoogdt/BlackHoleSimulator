use crate::curved_space;
use curved_space::Metric;
use curved_space::SpaceObject;

use crate::python_interface;
use crate::ray_tracer;
use crate::black_body;

extern crate image;
use self::image::{ImageBuffer, RgbImage};
use std::path::Path;


/// setup1: launch parallel photons in the direction of the black hole an watch it bend


// generate photons on line (x,y,z) = (-distance,0,i*spacing ), i=0..number going in the +x direction
fn generate_parallel_photons<'a>( distance : f64, spacing: f64, number: i32, metric : &'a  curved_space::SchwarzschildMetric ) -> Vec< curved_space::SchwarzschildObject<'a>  > {
    let mut a :Vec< curved_space::SchwarzschildObject  > = Vec::new(); 
    for i in 0..number {
        
        a.push( metric.spawn_space_object_from_cartesian( [ 0.0,-distance, 0.0,(i as f64)*spacing,],
          [1.0,1.0,0.0,0.0,] , 0.0 ))
    }
    return a;
}


pub fn launch_parallel_photons(){

    let num_photons = 50;
    let mut results : Vec<Vec<[f64;8]>> = Vec::new();

    let metric = curved_space::SchwarzschildMetric{r_s:1.0};

    let mut photons  = generate_parallel_photons(5.0,0.2,num_photons , & metric);

    for (i,x) in photons.iter_mut().enumerate() {

        let mut v : Vec<[f64;8]> = Vec::new();

        // println!("{}",x);
      
        for _ in 0..2000 {
            
            v.push( x.get_cartesian_coordinates_and_momenta() ); //creates a copy to push
            x.take_step(0.01);

            if x.get_coordinates()[1] < 1.05 {
                break;
            }
        
        } 

        let _ = python_interface::save_to_csv( &v, format!("files/photon{}.csv", i));

        results.push(v);
    }

    python_interface::launch_python(num_photons, "xz" ); 
}

//////////////////
/// 

pub fn ray_trace_schwarzshild(){
    let r_s = 1.0;
    let metric = curved_space::SchwarzschildMetric{ r_s : r_s };
    let camera = ray_tracer::Camera{ 
        pos : [0.0, -5.0,0.0,0.01],
        direction : [1.0,1.0,0.0,0.0],
        x_res : 800,
        y_res : 800,
        distance: 0.2,
        height : 0.3,
        width : 0.3, };

    let colored_sphere = ray_tracer::CollsionObject{ detect_collision : Box::new( 
            move |_: &[f64;4], new:  &[f64;4]|-> Option< image::Rgb<u8> > {
                if new[1] <= 1.05* r_s {
                    let p =  (10.0* new[3] / ( std::f64::consts::PI).round()  )as i64; 
                    let h =  (10.0 * new[2]  / ( std::f64::consts::PI).round()  )as i64; 
                    //println!("p{} h{} orig {} orig{}",p,h,new[2],new[3]);
                    if (p+h)%2  == 0 {
                        return Some( image::Rgb([255,0,0]) )
                    }
                    return  Some( image::Rgb([0,0,255]));
                }
                None
            }
        )};
        
    let outer_sphere = ray_tracer::CollsionObject{ detect_collision : Box::new( 
        move |_: &[f64;4], new:  &[f64;4]|-> Option< image::Rgb<u8> > {
            if new[1] >= 6.0* r_s {
                return  Some( image::Rgb([0,0,0]));     //todo: map picture
            }
            None
        }
    )};    

    let col_objects : Vec< ray_tracer::CollsionObject > = vec![colored_sphere, outer_sphere];

    let mut ray_tracer = ray_tracer::new( camera, &col_objects,  &metric,2000 ,false);

    
    ray_tracer.run_simulation(5);


    ray_tracer.generate_image("schw800x800.bmp");

    //ray_tracer.plot_paths();

}


pub fn ray_trace_minkowski(){
    let r_s = 1.0;

    let metric = curved_space::MinkowskiMetric{};
    let camera = ray_tracer::Camera{ 
        pos : [0.0, -5.0,0.0,0.01],
        direction : [1.0,1.0,0.0,0.0],
        x_res : 800,
        y_res : 800,
        distance: 0.2,
        height : 0.3,
        width : 0.3, };

    let colored_sphere = ray_tracer::CollsionObject{ detect_collision : Box::new( 
            move |_: &[f64;4], new:  &[f64;4]|-> Option< image::Rgb<u8> > {
                let r = ( new[1].powi(2) + new[2].powi(2)+new[3].powi(2)).sqrt();
                if r <= 1.05* r_s {

                    let theta = (new[3]/r).acos();
                    let phi =  new[2].atan2(new[1]);

                    let p =  (10.0*theta / ( std::f64::consts::PI).round()  )as i64; 
                    let h =  (10.0 * phi / ( std::f64::consts::PI).round()  )as i64; 

                    if (p + h)%2  == 0 {
                        return Some( image::Rgb([255,0,0]) )
                    }
                    return  Some( image::Rgb([0,0,255]));
                }
                None
            }
        )};
        
    let outer_sphere = ray_tracer::CollsionObject{ detect_collision : Box::new( 
        move |_: &[f64;4], new:  &[f64;4]|-> Option< image::Rgb<u8> > {
            let r = ( new[1].powi(2) + new[2].powi(2)+new[3].powi(2)).sqrt();
            if r >= 6.0* r_s {
                return  Some( image::Rgb([0,0,0]));     //todo: map picture
            }
            None
        }
    )};    

    let col_objects : Vec< ray_tracer::CollsionObject > = vec![colored_sphere, outer_sphere];

    let mut ray_tracer = ray_tracer::new( camera, &col_objects,  &metric,2000 ,false);

    
    ray_tracer.run_simulation(5);


    ray_tracer.generate_image("flat800x800.bmp");

    //ray_tracer.plot_paths();

}


// blackbodytester

pub fn generate_blackbody(T_min: f64, T_max:f64, steps: u32) {
    let cie = black_body::cie_lookup::new();
    

    let ysteps = ( (steps as f64) /10.0).ceil() as u32;

    let mut img: RgbImage = ImageBuffer::new(steps, ysteps);

    for i in 0..steps {
        let T = T_min + (T_max-T_min)*(i as f64)/(steps as f64);

        let col = cie.get_black_body_intensity(T, -200.0);
        for j in 0..ysteps {
            img.put_pixel(i, j , col );
        } 
        
    }

    let p = Path::new("blackbody.bmp");

    img.save_with_format( p,  image::ImageFormat::Bmp).unwrap();

}

fn test_wavelentgh_convo (){
    let cie = black_body::cie_lookup::new();

    let arr : [f64;3] = [435.8,546.1,700.0];

    for (i,x) in arr.iter().enumerate() {
        let col = cie.wavelength_to_rgb(*x);
        println!("lambda: {} rgb: {},{},{}",x ,col[0], col[1],col[2] );
    }

     
 }