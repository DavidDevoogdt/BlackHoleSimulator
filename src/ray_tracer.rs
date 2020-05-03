
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
   
    pub fn run_simulation(& mut self, max_tolerance: f64){
        self.photons.par_iter_mut().progress_count( (self.camera.x_res*self.camera.y_res) as u64 ).for_each(|photon| {
        //self.photons.iter_mut().for_each(|photon| {
            //println!("hello from {}", photon.steps_left );

            'outer: while photon.steps_left > 0 {
                photon.steps_left -=1;

                let d_lambda = photon.dynamic.get_d_lambda();
                photon.dynamic.take_step( d_lambda );

                let new_err = photon.dynamic.get_error_estimate();

                if new_err < max_tolerance {

                    photon.dynamic.sanitize_coordinates();
                    photon.dynamic.store_state();

                    let [coord,mom] = photon.dynamic.calculate_coordinates_and_contravariant_momenta();

                    for obj in photon.collision_object.iter() {
                        match  obj.detect_collision(&photon.prev_position,&coord,&mom   ) {
                            Some(x) => {
                                //photon.dynamic.print();
                                //println!("photon colided");
                                photon.final_color = Some(x); //todo account for redshift
                                photon.steps_left=0;
                                break 'outer;
                            },
                            None => { },
                        }
                    }

                    photon.prev_position = coord;

                    

                    if photon.save_path {
                    photon.path.push(  photon.dynamic.calculate_cartesian_coordinates_and_momenta() );  
                    }
                }else{
                    
                    photon.dynamic.restore_state();
                    photon.dynamic.set_d_lambda(Some( 0.5*photon.dynamic.get_d_lambda())  );
                }      
            }
        });

    }
    

    pub fn plot_paths( &self, mode: &str) {
        if self.save_path {

            for (i,photon) in self.photons.iter().enumerate() {

                //println!("{}",photon.path[0] [1]  );                
                let _ = python_interface::save_to_csv( &photon.path, format!("files/photon{}.csv", i)  );
            }

            let num = self.photons.len() as i32;

            println!("plotting {} photons in python", num);

            python_interface::launch_python( num , mode)

        }
    }

    pub fn generate_image(&self,name: &str) {
        let mut img: RgbImage = ImageBuffer::new(self.camera.x_res as u32, self.camera.y_res as u32);
        for x in 0..(self.camera.x_res as u32) {
            for y in 0..(self.camera.y_res as u32){
                //println!(" x:{} y:{}",x,y);

                let photon = &self.photons[ (y+x*(self.camera.y_res as u32))  as usize];

                let color =  match photon.final_color {
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
pub fn new<'a,  T: curved_space::Metric<'a> > ( camera: Camera,objects:  &'a Vec< Box< dyn CollsionObject> >,metric : &'a T, max_steps : i32 ,save_path : bool) -> RayTracer<'a,T> {
    let r =  (camera.direction[0].powi(2)+camera.direction[1].powi(2)+camera.direction[2].powi(2)).sqrt();

    let theta = (camera.direction[2]/r).acos();
    let phi =  camera.direction[1].atan2(camera.direction[0]);

    let dheigth = camera.height/ (camera.y_res as f64);
    let dwidth = camera.width/ (camera.x_res as f64);

    let mut alpha_x = theta.cos()*phi.cos()    ;
    let mut alpha_y = theta.cos()*phi.sin();
    let mut alpha_z = -theta.sin();

    let mut  beta_x = -phi.sin();
    let mut  beta_y = phi.cos();

    let alpha_norm =  (alpha_x .powi(2)+alpha_y .powi(2)+alpha_z .powi(2)).sqrt();
    let beta_norm =  (beta_x .powi(2)+beta_y .powi(2)).sqrt();

    alpha_x /= alpha_norm;
    alpha_y /= alpha_norm;
    alpha_z /= alpha_norm;

    beta_x /= beta_norm ;
    beta_y /= beta_norm ;

    let theta_x = (alpha_x*camera.rotation_angle.cos() - beta_x*camera.rotation_angle.sin())*(dheigth);
    let theta_y = (alpha_y*camera.rotation_angle.cos() - beta_y*camera.rotation_angle.sin())*(dheigth);
    let theta_z = (alpha_z*camera.rotation_angle.cos())*(dheigth);

    let phi_x = (alpha_x*(camera.rotation_angle.sin()) + beta_x*camera.rotation_angle.cos())*(dwidth);
    let phi_y = (alpha_y*(camera.rotation_angle.sin()) + beta_y*camera.rotation_angle.cos())*(dwidth);
    let phi_z = (alpha_z*(camera.rotation_angle.sin()))*(dwidth);

    //let r2 =  (camera.direction[1].powi(2) + camera.direction[2].powi(2) ).sqrt();


    let x_dir = camera.direction[0]*camera.distance/r;
    let y_dir = camera.direction[1]*camera.distance/r;
    let z_dir = camera.direction[2]*camera.distance/r;


    let mut photon_array = Vec::with_capacity( (camera.x_res*camera.y_res) as usize );

    for x in 0..camera.x_res{
        
        for y in 0..camera.y_res{

            let th =  ((2*y -camera.y_res + 1 ) as f64 ) /2.0 ;
            let ph = ((2*x -camera.x_res + 1 ) as f64 ) /2.0;


            let px = x_dir + th* theta_x+ ph*phi_x;
            let py = y_dir + th* theta_y+ ph*phi_y;
            let pz = z_dir + th* theta_z+ ph*phi_z;

            let norm =  (px.powi(2)+py.powi(2)+pz.powi(2)).sqrt();

            let pos = [0.0, camera.pos[0],camera.pos[1],camera.pos[2] ];

            let obj  = metric.spawn_space_object_from_cartesian( pos , [1.0, px/norm,  py/norm, pz/norm], 0.0 );
            

            let p = Photon{ 
                save_path: save_path,
                collision_object: objects ,
                prev_position : pos,
                dynamic: obj,
                steps_left: max_steps,
                path: vec![] ,
                phantom  :  std::marker::PhantomData,
                final_color: None
            };

            photon_array.push( p) ;
        }
    }
    
    RayTracer{ photons : photon_array, camera: camera,  save_path : save_path }   
}



//This creates the initial photons for the chosen metric and associated space object
pub fn new_parallel<'a,  T: curved_space::Metric<'a> > (metric : &'a T, objects:  &'a Vec< Box< dyn CollsionObject> >, direction:[f64;3] ,center: [f64;3], inc_vector: [f64;3], number_photons: i32,  max_steps : i32 ) -> RayTracer<'a,T> {
    
    //not relevant
    let cam = Camera{
        pos : center,
        direction : direction,
        x_res : number_photons,
        y_res :  1,
        distance: 1.0,
        width: 1.0,
        height : 0.0,
        rotation_angle : 0.0,
    };

    let mut photon_array = Vec::with_capacity( number_photons as usize );

    for x in 0..number_photons{

        let mut pos = [0.0, center[0],center[1],center[2] ];

        let inc_number = (-number_photons + 2*x )/2;

        for i in 0..3{
            pos[i+1] += (inc_number as f64)* inc_vector[i];
        }

        let norm =  (direction[0].powi(2)+direction[1].powi(2)+direction[2].powi(2)).sqrt();

        let obj  = metric.spawn_space_object_from_cartesian( pos , [1.0, direction[0]/norm,  direction[1]/norm, direction[2]/norm], 0.0 );
        
        let p = Photon{ 
            save_path: true,
            collision_object: objects,
            prev_position : pos,
            dynamic: obj,
            steps_left: max_steps,
            path: vec![] ,
            phantom  :  std::marker::PhantomData,
            final_color: None
        };

        photon_array.push( p) ;
        
    }
    
    RayTracer{ photons : photon_array, camera: cam,  save_path : true }   
}




//in cartesian coordinates
pub struct Camera {
    pub pos : [f64;3],
    pub direction : [f64;3],
    pub x_res : i32,
    pub y_res : i32,
    pub distance: f64,
    pub height : f64,
    pub width : f64,
    pub rotation_angle: f64,
}

pub trait CollsionObject: Send + Sync{
    fn  detect_collision (& self, old: &[f64;4], new: &[f64;4], momentum:  &[f64;4] ) -> Option < image::Rgb<u8>>;
}

//
pub struct Sphere {
    pub radius : f64,
    pub divisions : f64,
    pub color1 : image::Rgb<u8>,
    pub color2 : image::Rgb<u8>,
}

pub struct SphereCart {
    pub radius : f64,
    pub divisions : f64,
    pub color1 : image::Rgb<u8>,
    pub color2 : image::Rgb<u8>,
}


impl CollsionObject for Sphere{

    fn  detect_collision (& self, _: &[f64;4], new:  &[f64;4] , _:  &[f64;4] ) -> Option < image::Rgb<u8>>{  


        if new[1] <= self.radius {
            let p =  ( self.divisions * new[3] / ( std::f64::consts::PI).round()  )as i64; 
            let h =  ( self.divisions * new[2]  / ( std::f64::consts::PI).round()  )as i64; 
            //println!("p{} h{} orig {} orig{}",p,h,new[2],new[3]);
            if (p+h)%2  == 0 {
                return Some( self.color1 )
            }
            return  Some( self.color2 );
        }
        None
    }
}

impl CollsionObject for SphereCart{

    fn  detect_collision (& self, _: &[f64;4], new: &[f64;4], _:  &[f64;4]) -> Option < image::Rgb<u8>>{  


        let r = ( new[1].powi(2) + new[2].powi(2)+new[3].powi(2)).sqrt();
        if r <= self.radius {

            let theta = (new[3]/r).acos();
            let phi =  new[2].atan2(new[1]);

            let p =  (self.divisions *theta / ( std::f64::consts::PI).round()  )as i64; 
            let h =  (self.divisions * phi / ( std::f64::consts::PI).round()  )as i64; 

            if (p+h)%2  == 0 {
                return Some( self.color1 )
            }
            return  Some( self.color2 );
        }
        None
    }

}


// accretion disk

pub struct Annulus {
    pub radius1 : f64,
    pub radius2 : f64,
    pub divisions_radial : i32,
    pub divisions_angular : i32,
    pub color1 : image::Rgb<u8>,
    pub color2 : image::Rgb<u8>,
}

pub struct AnnulusCart {
    pub radius1 : f64,
    pub radius2 : f64,
    pub divisions_radial : i32,
    pub divisions_angular : i32,
    pub color1 : image::Rgb<u8>,
    pub color2 : image::Rgb<u8>,
}

impl CollsionObject for Annulus{

    fn  detect_collision (& self, old: &[f64;4], new: &[f64;4], _:  &[f64;4]) -> Option < image::Rgb<u8>>{  
        let dtheta1= std::f64::consts::PI / 2.0 - new[3];
        let dtheta2 = std::f64::consts::PI / 2.0 - old[3]; 



        //theta switches side
        if  dtheta1*dtheta2< 0.0 {
            if  (dtheta1-dtheta2).abs() <= 0.5 {
                let perc = dtheta2.abs()/ (dtheta1.abs()+dtheta2.abs());

                let r = perc*new[1] + (1.0-perc)*old[1];

                if r >= self.radius1 &&  r <= self.radius2 {  

                    let phi = perc*new[2] + (1.0-perc)*old[2]  +1.0;  
                   
                    let rfact  = ((self.divisions_radial as f64) * (r-self.radius1) / ( self.radius2-self.radius1))  as i32;
                    let p =  ((self.divisions_angular as f64) * phi / ( 2.0*std::f64::consts::PI))  as i32; 

                    if (p+rfact).rem_euclid(2)==0 {
                        return  Some( self.color1);
                    }

                    return  Some( self.color2);
                        
                }
            }           
        }
        None
    }
}

impl CollsionObject for AnnulusCart{

    fn  detect_collision (& self, old: &[f64;4], new: &[f64;4], _:  &[f64;4]) -> Option < image::Rgb<u8>>{  


        let r = ( new[1].powi(2) + new[2].powi(2)+new[3].powi(2)).sqrt();

        if r>= self.radius1 &&  r  <= self.radius2 {  

            let theta1 = (new[3]/r).acos();
            let theta2 = (old[3]/r).acos();
                

            //theta switches side
            if  (std::f64::consts::PI / 2.0 - theta1)* (std::f64::consts::PI / 2.0 - theta2) < 0.0 {
                
                let phi =  new[2].atan2(new[1]);

                let p =  (self.divisions_radial as f64 * phi / ( std::f64::consts::PI).round()  )as i64; 

                if p%2==0 {
                    return  Some( self.color1);
                }

                return  Some( self.color2);
                        
            }                
        }
        None
    }

}

// skybox

pub struct Skybox {
    pub x_res:i32,
    pub y_res:i32,
    pub image: image::RgbImage,
    pub radius: f64,
    pub phi_offset: f64,
} 

pub struct SkyboxCart {
    pub x_res:i32,
    pub y_res:i32,
    pub image: image::RgbImage,
    pub radius: f64,
    pub phi_offset: f64,
} 

pub struct SkyboxKerr {
    pub x_res:i32,
    pub y_res:i32,
    pub image: image::RgbImage,
    pub radius: f64,
    pub phi_offset: f64,
    pub a : f64,
} 

// equirectangular projection
impl CollsionObject for SkyboxKerr {

    fn  detect_collision (& self, _: &[f64;4], new: &[f64;4], momentum:  &[f64;4]) -> Option < image::Rgb<u8>>{  

        let r =new[1];

        if r > self.radius {

            // first needs to be covariant
            let [_,cart_mom] = curved_space::kerr_to_cart( &[new[1],new[2],new[3]] , & [momentum[1],momentum[2],momentum[3]], self.a );

            let mom_r = (cart_mom[0].powi(2) + cart_mom[1].powi(2)+cart_mom[2].powi(2)).sqrt();

            let theta = (cart_mom[2]/mom_r).acos();
            let phi = cart_mom[1].atan2( cart_mom[0] );


            let xcoor = (( (self.x_res as f64) * (  ( phi+self.phi_offset + std::f64::consts::PI )/(2.0* std::f64::consts::PI))).round() as i32  ).rem_euclid( self.x_res) ;
            let ycoor = (( (self.y_res as f64) * (  ( theta                )/(     std::f64::consts::PI))).round() as i32).rem_euclid( self.y_res);

            return Some( *self.image.get_pixel(xcoor as u32, ycoor as u32) )
            //return Some( image::Rgb([ (xcoor*255/self.x_res).try_into().unwrap() ,0,0]) )
        }

        None
    }
}


impl CollsionObject for Skybox {

    fn  detect_collision (& self, _: &[f64;4], new: &[f64;4], momentum:  &[f64;4]) -> Option < image::Rgb<u8>>{  

       
        let r =new[1];

        if r > self.radius {

            let [_,cart_mom] = curved_space::spher_to_cart( &[new[1],new[2],new[3]] , & [momentum[1],momentum[2],momentum[3]]);

            let mom_r = (cart_mom[0].powi(2) + cart_mom[1].powi(2)+cart_mom[2].powi(2)).sqrt();

            let theta = (cart_mom[2]/mom_r).acos();
            let phi = cart_mom[1].atan2( cart_mom[0] );

            let xcoor = (( (self.x_res as f64) * (  ( phi+self.phi_offset)/(2.0* std::f64::consts::PI))) as i32  ).rem_euclid( self.x_res) ;
            let ycoor = (( (self.y_res as f64) * (  ( theta                )/(     std::f64::consts::PI))) as i32 ) .rem_euclid( self.y_res);

            return Some( *self.image.get_pixel(xcoor as u32, ycoor as u32) )
            //return Some( image::Rgb([ (xcoor*255/self.x_res).try_into().unwrap() ,0,0]) )
        }

        None
    }
}

impl CollsionObject for SkyboxCart {
    fn  detect_collision (& self, _: &[f64;4], new: &[f64;4], momentum:  &[f64;4]) -> Option < image::Rgb<u8>>{  
       

        let r = ( new[1].powi(2) + new[2].powi(2)+new[3].powi(2)).sqrt();

        if r > self.radius{

            let cart_mom = momentum;

            let mom_r = (cart_mom[1].powi(2) + cart_mom[2].powi(2)+cart_mom[3].powi(2)).sqrt();

            let theta = (cart_mom[3]/mom_r).acos();
            let phi = cart_mom[2].atan2( cart_mom[1] );

            let xcoor = (( (self.x_res as f64) * (  ( phi+self.phi_offset)/(2.0* std::f64::consts::PI))) as i32  ).rem_euclid( self.x_res) ;
            let ycoor = (( (self.y_res as f64) * (  ( theta                )/(     std::f64::consts::PI))) as i32 ) .rem_euclid( self.y_res);

            return Some( *self.image.get_pixel(xcoor as u32, ycoor as u32) )
            //return Some( image::Rgb([ (xcoor*255/self.x_res).try_into().unwrap() ,0,0]) )
        }
        None
    }
}



//
pub struct Photon< 'a, T: curved_space::SpaceObject<'a> + std::marker::Send + std::marker::Sync >  {
    dynamic : T,
    prev_position :  [f64;4],
    steps_left: i32,
    phantom : std::marker::PhantomData<&'a T>,
    final_color : Option< image::Rgb<u8> >,
    collision_object : &'a Vec< Box< dyn CollsionObject> >,
    path : Vec< [[f64;4];2]>,
    save_path : bool,
}

