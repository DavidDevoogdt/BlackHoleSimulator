extern crate csv;

extern crate serde_derive;


use std::error::Error;
use csv::Writer;
use std::process::Command;

pub mod curved_space;

use curved_space::SpaceObject;
//use curved_space::SchwarzschildObject;



// generate photons on line (x,y,z) = (-distance,0,i*spacing ), i=0..number going in the +x direction
fn generate_parallel_photons<'a>( distance : f64, spacing: f64, number: i32, metric : &'a  curved_space::SchwarzschildMetric ) -> Vec< curved_space::SchwarzschildObject<'a>  > {
    let mut a :Vec< curved_space::SchwarzschildObject  > = Vec::new(); 
    for i in 0..number {
        a.push( curved_space::spawn_space_object_from_cartesian( [ 0.0,-distance, 0.0,(i as f64)*spacing,],
          [1.0,1.0,0.0,0.0,] , 0.0,  & metric ))
    }
    return a;
}

fn launch_python(n :i32){
    Command::new("python3")
            .arg("src/plotter.py")
            .arg(format!("-N {}",n) )
            .spawn()
            .expect("error while calling python");

}

fn save_to_csv(v : &Vec< [f64;8] >, name : String) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(  name )?;

    for x in v {
        wtr.serialize( x )?;
    }

    wtr.flush()?;
    Ok(())
}


fn main() {



    let num_photons = 50;
    let mut results : Vec<Vec<[f64;8]>> = Vec::new();

    let metric = curved_space::SchwarzschildMetric{r_s:1.0};

    let mut photons  = generate_parallel_photons(5.0,0.2,num_photons , & metric);

    for (i,x) in photons.iter_mut().enumerate() {

        let mut v : Vec<[f64;8]> = Vec::new();

        // println!("{}",x);
      
        for _ in 0..2000 {
            
            v.push( x.get_cartesian_coordinates_and_momenta() ); //creates a copy to push
            x.take_step();

            if x.get_coordinates()[1] < 1.05 {
                break;
            }
        
        } 

        let _ = save_to_csv( &v, format!("files/photon{}.csv", i));

        results.push(v);
    }

     launch_python(num_photons); 
 }

