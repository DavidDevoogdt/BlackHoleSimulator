use crate::curved_space;
use curved_space::Metric;
use curved_space::SpaceObject;

use crate::python_interface;
use crate::ray_tracer;

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

    python_interface::launch_python(num_photons); 
}

//////////////////
/// 

pub fn ray_trace(){
    let metric = curved_space::SchwarzschildMetric{ r_s : 1.0 };
    let camera = ray_tracer::Camera{ 
        pos : [0.0, -5.0,0.0,0.0],
        direction : [1.0,1.0,0.0,0.0],
        x_res : 8,
        y_res : 8,
        distance: 0.5,
        height : 0.5,
        width : 0.5, };
    let col_objects : Vec< ray_tracer::CollsionObject > = vec![];

    let mut ray_tracer = ray_tracer::new( camera, col_objects,  &metric);

    for _ in 0..10 {
        ray_tracer.take_step();
    }
}

