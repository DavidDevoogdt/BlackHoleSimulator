
extern crate image;
use self::image::{ImageBuffer, RgbImage};
use std::path::Path;

extern crate rayon;



use crate::curved_space;
use curved_space::SpaceObject;
use crate::python_interface;

extern crate indicatif;
use self::indicatif::ParallelProgressIterator;
use self::rayon::iter::{ParallelIterator, IntoParallelRefMutIterator };


pub struct RayTracer<'a,T: curved_space::Metric<'a>  > {
    camera : Camera,
    photons : Vec<Photon<'a, T::SpecificSpaceObject> > ,
    save_path : bool,
}


impl<'a,T: curved_space::Metric<'a>  > RayTracer<'a,T> {
   
    pub fn run_simulation(& mut self, steps_per_collision_detection: u8){
        self.photons.par_iter_mut().progress_count( (self.camera.x_res*self.camera.y_res) as u64 ).for_each(|photon| {
            //println!("hello from {}", photon.steps_left );
            let mut steps_taken : u8 = 0; 
            'outer: while photon.steps_left > 0 {
                 photon.steps_left -=1;

                 photon.dynamic.take_step( photon.d_lambda );
                 steps_taken +=1;

                if steps_taken == steps_per_collision_detection {
                     for obj in photon.collision_object.iter() {
                        match (*obj.detect_collision )( &photon.prev_position , photon.dynamic.get_coordinates() ) {
                            Some(x) => {
                                //println!("photon colided");
                                photon.final_color = Some(x); //todo account for redshift
                                photon.steps_left=0;
                                break 'outer;
                            },
                            None => { },
                        }
                     }

                     steps_taken = 0;
                 }

                if photon.save_path {
                   photon.path.push(  photon.dynamic.get_cartesian_coordinates_and_momenta() );  
                }
            }
        });

    }
    

    pub fn plot_paths( &self) {
        if self.save_path {

            for (i,photon) in self.photons.iter().enumerate() {
                let _ = python_interface::save_to_csv( &photon.path, format!("files/photon{}.csv", i)  );
            }

            let num = self.photons.len() as i32;

            println!("plotting {} photons in python", num);

            python_interface::launch_python( num , "xyz")

        }
    }

    pub fn generate_image(&self,name: &str) {
        let mut img: RgbImage = ImageBuffer::new(self.camera.x_res as u32, self.camera.y_res as u32);
        for x in 0..(self.camera.x_res as u32) {
            for y in 0..(self.camera.y_res as u32){
                //println!(" x:{} y:{}",x,y);
                let color =  match self.photons[ (x+y*(self.camera.x_res as u32))  as usize].final_color {
                    None => image::Rgb([0,255,0]),
                    Some(x) =>x, 
                } ;

                img.put_pixel(x, y , color  );
            }
        }


        let p = Path::new(name);

        img.save_with_format( p,  image::ImageFormat::Bmp).unwrap();
        
    }
}


//This creates the initial photons for the chosen metric and associated space object
pub fn new<'a,  T: curved_space::Metric<'a> > ( camera: Camera, objects : &'a Vec<CollsionObject>, metric : &'a T, max_steps : i32 ,save_path : bool ) -> RayTracer<'a,T> {
    let r = ( camera.direction[1].powi(2) + camera.direction[2].powi(2) + camera.direction[3].powi(2)).sqrt();

    let theta = (camera.direction[3]/r).acos();
    let phi =  camera.direction[2].atan2(camera.direction[1]);

    //let r2 =  (camera.direction[1].powi(2) + camera.direction[2].powi(2) ).sqrt();


    let dheigth = camera.height/ (camera.x_res as f64);
    let dwidth = camera.width/ (camera.y_res as f64);

    let dtheta = ( dheigth / camera.distance).atan();
    let dphi = ( dwidth / camera.distance).atan();

    let mut photon_array = Vec::with_capacity( (camera.x_res*camera.y_res) as usize );

    for x in 0..camera.x_res{
        
        for y in 0..camera.y_res{


            //println!( "{}", ((2*y -camera.y_res + 1 ) as f64 ) /2.0 );
 
            let th = theta + ((2*y -camera.y_res + 1 ) as f64 ) /2.0 *dtheta;
            let ph = phi+ ((2*x -camera.x_res + 1 ) as f64 ) /2.0 *dphi;
            let thc = th.cos();
            let ths = th.sin();
            let phc = ph.cos();
            let phs = ph.sin();

            let obj  = metric.spawn_space_object_from_cartesian( camera.pos , [1.0, ths*phc, ths*phs, thc], 0.0 );
            
            let prev_pos = obj.get_coordinates().clone();

            let p = Photon{ save_path: save_path, collision_object: objects ,prev_position : prev_pos,dynamic: obj , d_lambda: 0.01, steps_left: max_steps,path: vec![] ,phamtom  :  std::marker::PhantomData, final_color: None };

            photon_array.push( p) ;
        }
    }
    
    RayTracer{ photons : photon_array, camera: camera,  save_path : save_path }   
}


//in cartesian coordinates
pub struct Camera {
    pub pos : [f64;4],
    pub direction : [f64;4],
    pub x_res : i32,
    pub y_res : i32,
    pub distance: f64,
    pub height : f64,
    pub width : f64,
}

pub struct CollsionObject { //collision
    // first old coor, then new ones
    pub detect_collision :  Box< dyn Fn( &[f64;4], &[f64;4]) -> Option < image::Rgb<u8>>  + Send + Sync >  ,  //gets coordinates as input and returns option with color in wavelenth 
    //bounding_box : [f64;3],
}

pub struct Photon< 'a, T: curved_space::SpaceObject<'a> + std::marker::Send + std::marker::Sync >  {
    dynamic : T,

    prev_position :  [f64;4],
    d_lambda : f64,
    steps_left: i32,
    phamtom : std::marker::PhantomData<&'a T>,
    final_color : Option< image::Rgb<u8> >,
    collision_object : &'a Vec< CollsionObject>,
    path : Vec< [f64;8]>,
    save_path : bool,
}

