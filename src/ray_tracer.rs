
use crate::curved_space;
use curved_space::SpaceObject;

pub struct RayTracer<'a,T: curved_space::Metric<'a> > {
    metric: &'a T,
    camera : Camera,
    objects : Vec<CollsionObject>,
    photons : Vec<Photon<'a, T::SpecificSpaceObject> > ,

}


impl<'a,T: curved_space::Metric<'a>> RayTracer<'a,T> {
    pub fn take_step(&mut self){
        for  photon in self.photons.iter_mut() {
            photon.dynamic.take_step( photon.d_lambda );
            photon.steps_left -=1;
        }
    }   
}


//This creates the initial photons for the chosen metric and associated space object
pub fn new<'a,  T: curved_space::Metric<'a> > ( camera: Camera, objects : Vec<CollsionObject>, metric : &'a T ) -> RayTracer<'a,T> {
    let r = ( camera.direction[1].powi(2) + camera.direction[2].powi(2) + camera.direction[3].powi(2)).sqrt();

    let theta = (camera.direction[3]/r).acos();
    let phi =  camera.direction[2].atan2(camera.direction[1]);

    //let r2 =  (camera.direction[1].powi(2) + camera.direction[2].powi(2) ).sqrt();


    let dheigth = camera.height/ (camera.x_res as f64);
    let dwidth = camera.width/ (camera.y_res as f64);

    let dtheta = ( dheigth / camera.distance).atan();
    let dphi = ( dwidth / camera.distance).atan();

    let mut photon_array = Vec::with_capacity( (camera.x_res*camera.y_res) as usize );

    for x in -camera.x_res/2+1..camera.x_res/2{
        
        for y in -camera.y_res/2+1..camera.y_res/2{
            let th = theta + (y as f64)*dtheta;
            let ph = phi+(x as f64)*dphi;
            let thc = th.cos();
            let ths = th.sin();
            let phc = ph.cos();
            let phs = ph.sin();

            let obj  = metric.spawn_space_object_from_cartesian( camera.pos , [1.0, ths*phc, ths*phs, thc], 0.0 );
            let mom_copy = obj.get_momenta().clone();

            let p = Photon{ dynamic: obj , d_lambda: 0.01, steps_left: 200, initial_momentum: mom_copy, phamtom  :  std::marker::PhantomData };

            photon_array[ (x+y*camera.x_res) as usize] = p;
        }
    }
    
    RayTracer{ photons : photon_array, camera: camera, metric: metric, objects: objects }   
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
    detect_collision : Box< dyn Fn( &[f64;4], &[f64;4],f64 ) -> Option<f64> >,  //gets coordinates as input and returns option with color in wavelenth 
    bounding_box : [f64;3],
}

pub struct Photon< 'a, T: curved_space::SpaceObject<'a>>  {
    dynamic : T,
    initial_momentum : [f64;4], //used to calculate the energy and place on the screen
    d_lambda : f64,
    steps_left: i32,
    phamtom : std::marker::PhantomData<&'a T>,
}

