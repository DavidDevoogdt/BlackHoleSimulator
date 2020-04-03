// Helper objects 

/// integrate equation dy/dt = f(t,y) over a time step dt and store result in orig params
pub fn rk4 (f : impl Fn(f64) -> f64, y: &mut f64,dt:f64)  {
    let k1 = dt* f(*y);
    let k2 = dt* f(*y+k1/2.0);
    let k3 = dt* f(*y+k2/2.0);
    let k4 = dt* f(*y+k3);
    *y=*y+(k1+2.0*k2+2.0*k3+k4)/6.0;
}


#[test]
fn integrator(){
    let mut t= 0.0;
    let mut y= 1.0;
    for _ in 0..100 {
        let dt = 0.01;
         rk4( |y: f64| -> f64 {y} ,&mut y,dt)  ;
        t = t+ dt;
        println!("t:{:.3} y:{:.3}  e^t:{:.3}",t,y,t.exp() );
    }
}


///// main definitions
pub trait SpaceObject{
    fn get_metric(&self) -> & dyn Metric;
    fn get_coordinates(&self) -> &[f64;4];
    fn get_momenta(&self) -> &[f64;4];

    fn get_mut_coordinates(&mut self) -> &mut [f64;4];
    fn get_mut_momenta(&mut self) -> &mut [f64;4];

    fn get_coordinates_and_momenta(&self) -> [f64;8] {
        let coor = self.get_coordinates();
        let mom = self.get_momenta(); 
        [coor[0], coor[1],coor[2], coor[3],mom[0], mom[1],mom[2], mom[3], ]
    }
    //default implemented things
    fn get_contravariant_momenta (&self) -> [f64;4]{
        self.get_metric().to_contra( self.get_coordinates(), self.get_momenta() )    
    }

    // infer p^t from p^mu p_mu = (m*c)^2
    fn contract_momentum(&self ) -> f64{
        let p1 = self.get_momenta();
        let p2 = self.get_contravariant_momenta();

        let mut sum = 0.0;

        for i in 0..3 {
            sum += p1[i] + p2[i];
        }
        return sum;
    }


    fn take_step(&mut self){

        let copy_coor =   self.get_coordinates().clone();
        let copy_mom =   self.get_momenta().clone();

        //closures for the RK4 integration based on geodesic equation  (Dp^mu/D tau = 0 for force free field)
        let dp_r =  self.get_metric().covariant_derivative(1, &copy_coor, &copy_mom  ) ;
        let dp_h =  self.get_metric().covariant_derivative(3, &copy_coor, &copy_mom  ) ;
        let dp_t =  self.get_metric().covariant_derivative(0, &copy_coor,&copy_mom  ) ;
        let dp_p =  self.get_metric().covariant_derivative(2, &copy_coor,&copy_mom   ) ;

        let d_lambda = 0.01; //lambda is defined as tau/m

        let mut_mom = self.get_mut_momenta();

        //update the momenta 
        rk4( dp_r, &mut mut_mom[1] , d_lambda);
        rk4( dp_h, &mut mut_mom[3] , d_lambda);
        rk4( dp_p, &mut mut_mom[2] , d_lambda);
        rk4( dp_t, &mut mut_mom[0] , d_lambda);

        let mut_coord = self.get_mut_coordinates();
    
        //update coordinates
        rk4( |_| copy_mom[1], &mut mut_coord[1] , d_lambda);
        rk4( |_| copy_mom[2], &mut mut_coord[2] , d_lambda);
        rk4( |_| copy_mom[3], &mut mut_coord[3], d_lambda);
        rk4( |_| copy_mom[0], &mut mut_coord[0] , d_lambda);

        //todo check for collision with target
       

    }
}


pub trait Metric {

    fn get_christoffel (&self,x0: i8,x1: i8,x2: i8,coor: &[f64;4]) -> f64;
    fn g_lower(&self, x0 :i8, x1 :i8,coor :&[f64;4]) -> f64;

    fn to_contra(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4];
    fn covariant_derivative <'a>(&self,index:i8, coor: &'a [f64;4],  vec :&'a [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>; // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 

}

////////////////////implementation of the schwarzschildmetric

//t,r,phi,theta
pub struct SchwarzschildMetric {
    pub r_s : f64,
}

impl Metric for SchwarzschildMetric {

    fn get_christoffel (&self,x0: i8,x1: i8,x2: i8,coor: &[f64;4]) -> f64{
        match (x0,x1,x2){
            (0,1,0 ) => self.r_s/ (2.0*coor[1]*(coor[1]-self.r_s)),
            (0,0,0) => -self.r_s/ (2.0*coor[1]*(coor[1]-self.r_s)),
            (1,0,0) => self.r_s*(coor[1]-self.r_s)/(2.0* coor[1].powi(3)),
            (1,2,2) => -(coor[1]-self.r_s)*  coor[3].sin() .powi(2),
            (1,3,3) => -(coor[1]-self.r_s),
            (3,1,3) => 1.0/coor[1],
            (2,1,2) => 1.0/coor[1],
            (3,2,2) => -coor[3].sin()*coor[3].cos(),
            (2,3,2) => coor[3].cos()/coor[3].sin(),
            _ => 0.0,
        }
    }

    fn g_lower(&self, x0 :i8, x1 :i8,coor :&[f64;4]) -> f64{
        match (x0,x1) {
            (0,0) => -(1.0 - self.r_s/coor[ 1 ] ),
            (1,1) => 1.0/(1.0 - self.r_s/coor[ 1 ] ),
            (2,2) => coor[1].powi(2),
            (3,3) => (coor[1] * coor[3].sin()  ).powi(2),
            __ => 0.0,
        }
    }

    fn to_contra(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ self.g_lower( 0, 0, coordinates )* vec[0], self.g_lower( 1, 1, coordinates )* vec[1]  ,self.g_lower( 1, 1, coordinates )* vec[1],self.g_lower( 1, 1, coordinates )* vec[1]]
    }

    // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 
    fn covariant_derivative<'a> (&self,index:i8, coor: &'a [f64;4], vec :&'a [f64;4]  ) ->  Box<dyn Fn(f64)->f64+ 'a>{

        match index {
            0=>{let trt = self.get_christoffel(0,1,0, coor);
                Box::new( move |x_0: f64| -> f64 { -(2.0* trt *  vec[1]* x_0)} )},
            1=>{
                let rrr = self.get_christoffel(1,1,1, coor);
                let rtt = self.get_christoffel(1,0,0, coor);
                let rpp = self.get_christoffel(1,2,2, coor);
                let rhh = self.get_christoffel(1,3,3, coor);

                Box::new(move |x_1: f64| -> f64 { -(rrr* x_1.powi(2) +  rtt *  vec[0].powi(2) + rpp *  vec[2].powi(2) + rhh *  vec[3].powi(2)   )})
            },
            2=> {
                let prp = self.get_christoffel(2,1,2, coor);
                let php = self.get_christoffel(2,3,2, coor);

                Box::new(move |x_2: f64| -> f64 { -(2.0*prp* vec[1]* x_2 + 2.0*php* vec[3]* x_2)})
            },
            3=> {
                let hrh = self.get_christoffel(3,1,3, coor);
                let hpp = self.get_christoffel(3,2,2, coor);
                Box::new(move |x_3: f64| -> f64 { -(2.0*hrh* vec[1]* x_3 + hpp* vec[2].powi(2))})
            },
                _ => { Box::new( move |_: f64| -> f64 { 0.0})},
        }
    } 
}

//implementation of a object living in a schwarzschildmetric

pub struct SchwarzschildObject <'a> {
    coordinates : [f64;4],
    momenta :  [f64;4],
    metric: &'a SchwarzschildMetric,
    _mass : f64,
}

impl<'a> SpaceObject for SchwarzschildObject<'a>{

    fn get_metric(& self) -> & dyn Metric {
        self.metric as & dyn Metric         
    }
    fn get_coordinates(&self) -> &[f64;4]{
        &self.coordinates
    }
    fn get_momenta(&self) -> &[f64;4]{
        &self.momenta
    }
    fn get_mut_coordinates(&mut self) -> &mut [f64;4]{
        &mut self.coordinates
    }
    fn get_mut_momenta(&mut self) -> &mut [f64;4]{
        &mut self.momenta
    }
}

impl<'a> SchwarzschildObject<'a>{
    pub fn get_cartesian_coordinates_and_momenta(&self) -> [f64;8] {
        let coor = self.get_coordinates();
        let mom = self.get_momenta(); 
    
        let ps = coor[2].sin();
        let pc = coor[2].cos();
        let hs = coor[3].sin();
        let hc = coor[3].cos();
        
        let pr = mom[1];
        let ph = mom[3];
        let pp = mom[2];
        
        let r = coor[1];
        
        [coor[0], r* hs*pc,r* hs*ps,r* hc,
        mom[0], 
        hs*pc*pr + r*hc*pc*ph - r*hs*ps*pp,
        hs*ps*pr + r*hc*ps*ph + r*ps*pc*pp,
        hc*pr -r*hs*ph]
        //todo convert time and energy properly
    }
    
}


pub fn spawn_space_object<'a>( coordinates : [f64;4], momenta : [f64;4] , mass : f64, metric : &'a SchwarzschildMetric ) -> SchwarzschildObject<'a>{ 
    SchwarzschildObject{ coordinates : coordinates, momenta: momenta, _mass: mass, metric: metric  }
}

pub fn spawn_space_object_from_cartesian<'a>( coordinates :[f64;4], momenta : [f64;4] , mass : f64, metric : &'a SchwarzschildMetric ) -> SchwarzschildObject<'a>{

    let r = ( coordinates[1].powi(2) + coordinates[2].powi(2) + coordinates[3].powi(2)).sqrt();

    let theta = (coordinates[3]/r).acos();
    let phi =  coordinates[2].atan2(coordinates[1]);

    let r2 =  (coordinates[1].powi(2) + coordinates[2].powi(2) ).sqrt();

    let pc = coordinates[1]/r2;
    let ps = coordinates[2]/r2;
    let hc = coordinates[3]/r;
    let hs = r2 /r;

    let newpos = [
        coordinates[0],
        r,
        phi,
        theta];

    let newmom =  [ 
    momenta[0], //todo infer this momentum from energy mass relation
    hs*pc* momenta[1] + hs*ps*momenta[2]+ hc*momenta[3],
    -ps/(r*hs)*momenta[1] + pc/(r*hs) *momenta[2],
    hc*pc/r *momenta[1] + hc*ps/r*momenta[2]-hs/r *momenta[3]];


    spawn_space_object( newpos , newmom, mass, metric )

}


////////////////////implementation of minkowski metric 