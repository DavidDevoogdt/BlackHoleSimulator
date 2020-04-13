

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
pub trait SpaceObject<'a>{
    //fn get_metric(&self) -> &dyn Metric<SpecificSpaceObject = Self>;
    fn get_coordinates_patch(&self) -> &[f64;4];
    fn get_momenta_patch(&self) -> &[f64;4];



    fn get_mut_coordinates_patch(&mut self) -> &mut [f64;4];
    fn get_mut_momenta_patch(&mut self) -> &mut [f64;4];

    // for spherical coordinates, multiple patches could be useful
    fn sanitize_coordinates(&mut self);


    fn get_coordinates_and_momenta(&self) -> [[f64;4];2] ;


    //
    fn get_cartesian_coordinates_and_momenta(&self) -> [f64;8];

    // copy paste these implementations in the corresponding section 
    // reason: size of Self of metric is unknown for generic implementation due to associated type

    fn get_contravariant_momenta (&self) -> [f64;4];
    fn covariant_derivative (&self,index:i8, coor: & [f64;4],  vec :& [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>;

    /*fn get_contravariant_momenta (&self) -> [f64;4]{
        self.get_metric().to_contra( self.get_coordinates(), self.get_momenta() )    
    }

    fn covariant_derivative (&'a self,index:i8, coor: &'a [f64;4],  vec :&'a [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>{
        self.get_metric().covariant_derivative(index, coor,vec )
    }*/
   
    
    fn contract_momentum(&self ) -> f64{
        let p1 = self.get_momenta_patch();
        let p2 = self.get_contravariant_momenta();

        let mut sum = 0.0;

        for i in 0..3 {
            sum += p1[i] + p2[i];
        }
        return sum;
    }


    fn take_step(&mut self, d_lambda : f64){

        let copy_coor =   self.get_coordinates_patch().clone();
        let copy_mom =   self.get_momenta_patch().clone();

        //closures for the RK4 integration based on geodesic equation  (Dp^mu/D tau = 0 for force free field)
        let dp_r =  self.covariant_derivative(1, &copy_coor, &copy_mom  ) ;
        let dp_h =  self.covariant_derivative(3, &copy_coor, &copy_mom  ) ;
        let dp_t =  self.covariant_derivative(0, &copy_coor,&copy_mom  ) ;
        let dp_p =  self.covariant_derivative(2, &copy_coor,&copy_mom   ) ;

        let mut_mom = self.get_mut_momenta_patch();

        //update the momenta 
        rk4( dp_r, &mut mut_mom[1] , d_lambda);
        rk4( dp_h, &mut mut_mom[3] , d_lambda);
        rk4( dp_p, &mut mut_mom[2] , d_lambda);
        rk4( dp_t, &mut mut_mom[0] , d_lambda);

        let mut_coord = self.get_mut_coordinates_patch();
    
        //update coordinates
        rk4( |_| copy_mom[1], &mut mut_coord[1] , d_lambda);
        rk4( |_| copy_mom[2], &mut mut_coord[2] , d_lambda);
        rk4( |_| copy_mom[3], &mut mut_coord[3], d_lambda);
        rk4( |_| copy_mom[0], &mut mut_coord[0] , d_lambda);

        //todo check for collision with target
       

    }
}


pub trait Metric<'a> : std::marker::Sync + std::marker::Send {
    type SpecificSpaceObject : SpaceObject<'a> + std::marker::Sync + std::marker::Send ;

    fn get_christoffel (&self,x0: i8,x1: i8,x2: i8,coor: &[f64;4]) -> f64;
    fn g_lower(&self, x0 :i8, x1 :i8,coor :&[f64;4]) -> f64;

    fn to_contra(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4];
    fn covariant_derivative (&self,index:i8, coor: & [f64;4],  vec :& [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>; // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 

    fn spawn_space_object(&'a self, coordinates : [f64;4], momenta : [f64;4] , mass : f64 ) -> Self::SpecificSpaceObject;
    fn spawn_space_object_from_cartesian(&'a self, coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject;

}

////////////////////implementation of the schwarzschildmetric

//t,r,phi,theta
pub struct SchwarzschildMetric {
    pub r_s : f64,
}

impl<'a> Metric<'a> for SchwarzschildMetric {
    type SpecificSpaceObject= SchwarzschildObject<'a>;

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
        [ self.g_lower( 0, 0, coordinates )* vec[0], self.g_lower( 1, 1, coordinates )* vec[1]  ,self.g_lower( 2, 2, coordinates )* vec[2],self.g_lower( 3, 3, coordinates )* vec[3] ]
    }

    // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 
    fn covariant_derivative (&self,index:i8, coor: & [f64;4], vect :& [f64;4]  ) ->  Box<dyn Fn(f64)->f64+ 'a>{

        let vec = vect.clone();

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

    fn spawn_space_object(&'a self ,coordinates : [f64;4], momenta : [f64;4] , mass : f64 ) -> Self::SpecificSpaceObject{ 
        let mut  res = SchwarzschildObject{ patch_coordinates : coordinates, patch_momenta: momenta, _mass: mass, metric: &self, coordinate_patch: 0, theta_threshold: 0.6*std::f64::consts::PI/2.0 };
        
        res.sanitize_coordinates();

        res
    }



    fn spawn_space_object_from_cartesian( &'a self ,coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{

        //println!("spawning {}{}{}{}",coordinates[0],coordinates[1],coordinates[2],coordinates[3]);

        let [coor,mom] = cart_to_spher(   &[coordinates[1],coordinates[2],coordinates[3]] ,  &[momenta[1],momenta[2],momenta[3]] );

        //println!("new pos {};{};{};",coor[0],coor[1],coor[2]);

        let newpos = [
            coordinates[0],
            coor[0],
            coor[1],
            coor[2]
        ];

        let newmom =  [ 
            momenta[0], //todo infer this momentum from energy mass relation
            mom[0],
            mom[1],
            mom[2]
        ];


        self.spawn_space_object( newpos , newmom, mass)

    }
}

//implementation of a object living in a schwarzschildmetric

///idea; the main soherical coordinates are XYZ with theta angle between Z axis and vector
/// if angle is smaller than threshold, another coordinate patch is used to improve numerical accuracy
/// XYZ->YZX
/// this mode is used for simulation only and the get_coordinates ruteruns the correct version 
pub struct SchwarzschildObject <'a> {
    patch_coordinates : [f64;4],
    patch_momenta :  [f64;4],
    metric: &'a SchwarzschildMetric,
    _mass : f64,
    coordinate_patch: u8, 
    theta_threshold: f64,
}

impl<'a> SpaceObject<'a> for SchwarzschildObject<'a>{

    // fn get_metric(& self) -> & dyn Metric<SpecificSpaceObject = SchwarzschildObject > {
    //     self.metric as & dyn Metric<SpecificSpaceObject = SchwarzschildObject >       
    // }
    fn get_coordinates_patch(&self) -> &[f64;4]{
        &self.patch_coordinates
    }
    fn get_momenta_patch(&self) -> &[f64;4]{
        &self.patch_momenta
    }
    fn get_mut_coordinates_patch(&mut self) -> &mut [f64;4]{
        &mut self.patch_coordinates
    }
    fn get_mut_momenta_patch(&mut self) -> &mut [f64;4]{
        &mut self.patch_momenta
    }

    fn get_coordinates_and_momenta(&self) -> [[f64;4];2] {
        match self.coordinate_patch {
            0=> [*self.get_coordinates_patch(), *self.get_momenta_patch()]  ,
            1=> { 
                let res = self.rotate_coordinate_patch();
                [
                    [ self.patch_coordinates[0], res[0][0],res[0][1],res[0][2] ],
                    [self.patch_momenta[0], res[1][0],res[1][1],res[1][2]],
                ]
            },
            _=> panic!("invalid patch state")
        }

     
    }

  

    fn get_contravariant_momenta (&self) -> [f64;4]{
        self.metric.to_contra( self.get_coordinates_patch(), self.get_momenta_patch() )    
    }

    fn covariant_derivative (& self,index:i8, coor: & [f64;4],  vec :& [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>{
        self.metric.covariant_derivative(index, coor,vec )
    }

    fn get_cartesian_coordinates_and_momenta(&self) -> [f64;8] {
        let mut ret : [f64;8] = [0.0;8];
        let [coordinates,momenta] = self.get_coordinates_and_momenta();

        let [coor,mom] = spher_to_cart(  &[coordinates[1],coordinates[2],coordinates[3]] ,  &[momenta[1],momenta[2],momenta[3]]);

        

        ret[0] = coordinates[0];
        ret[4] = momenta[0];

        for i in 1..4 {
            ret[i] = coor[i-1];
        }

        for i in 5..8 {
            ret[i] = mom[i-5];
        }
        
        ret
    }

    fn sanitize_coordinates(&mut self) {
        self.patch_coordinates[2] = self.patch_coordinates[2]%(std::f64::consts::PI*2.0);

        //x,y,z -> yzx
        if  (std::f64::consts::PI/2.0 - self.patch_coordinates[3]).abs() > self.theta_threshold {
            
            self.switch_patch();
        }
    }
}

//XYZ->YZX

impl<'a> SchwarzschildObject<'a>{
    fn rotate_coordinate_patch(&self)->  [[f64;3];2]  {
        
        let mut_coor = self.get_coordinates_patch();
        let mut_mom = self.get_momenta_patch();

        let [coor_cart,mom_cart] = spher_to_cart( &[mut_coor[ 1  ],mut_coor[ 2 ],mut_coor[3]] , &[mut_mom[1],mut_mom[2],mut_mom[3]] );

        match self.coordinate_patch{
            0 => {
                cart_to_spher(  &[ coor_cart[ 2 ],coor_cart[0],coor_cart[1]], &[mom_cart[2],mom_cart[0],mom_cart[1]])
            },
            1 => {
                cart_to_spher(  &[ coor_cart[ 1 ],coor_cart[2],coor_cart[0]], &[mom_cart[1],mom_cart[2],mom_cart[0]])
            },
            _=>panic!("invalid coodinate patch")
        }  
    }

    pub fn switch_patch(&mut self){
        //print!("rotated!");
            
        let [coor,mom] = self.rotate_coordinate_patch();

        let mut_coor = self.get_mut_coordinates_patch();

        mut_coor[1] = coor[0];
        mut_coor[2] = coor[1];
        mut_coor[3] = coor[2];

        let mut_mom = self.get_mut_momenta_patch();

        mut_mom[1] =  mom[0];
        mut_mom[2] =  mom[1];
        mut_mom[3] =  mom[2];

        self.coordinate_patch = (self.coordinate_patch + 1)%2;
    }

}


////////////////////implementation of minkowski metric 




//t,x,y,z
pub struct MinkowskiMetric {
}

impl<'a> Metric<'a> for MinkowskiMetric {
    type SpecificSpaceObject= MinkowskiObject<'a>;

    fn get_christoffel (&self,_: i8,_: i8,_ :i8,_: &[f64;4]) -> f64{
        return  0.0;
    }

    fn g_lower(&self, x0 :i8, x1 :i8,_ :&[f64;4]) -> f64{
        match (x0,x1) {
            (0,0) => 1.0 ,
            (1,1) => -1.0,
            (2,2) => -1.0,
            (3,3) => -1.0,
            __ => 0.0,
        }
    }

    fn to_contra(&self,_ : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ vec[0], - vec[1]  ,- vec[2],- vec[3]]
    }

    // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 
    fn covariant_derivative (&self,_:i8, _: & [f64;4], _:& [f64;4]  ) ->  Box<dyn Fn(f64)->f64+ 'a>{
         Box::new( move |_: f64| -> f64 { 0.0})
    } 

    fn spawn_space_object(&'a self ,coordinates : [f64;4], momenta : [f64;4] , mass : f64 ) -> Self::SpecificSpaceObject{ 
        MinkowskiObject{ coordinates : coordinates, momenta: momenta, _mass: mass, metric: &self }
    }

    fn spawn_space_object_from_cartesian( &'a self ,coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{
        self.spawn_space_object( coordinates , momenta, mass)
    }
}

//implementation of a object living in a schwarzschildmetric

pub struct MinkowskiObject <'a> {
    coordinates : [f64;4],
    momenta :  [f64;4],
    metric: &'a MinkowskiMetric,
    _mass : f64,
}

impl<'a> SpaceObject<'a> for MinkowskiObject<'a>{

    fn get_coordinates_patch(&self) -> &[f64;4]{
        &self.coordinates
    }
    fn get_momenta_patch(&self) -> &[f64;4]{
        &self.momenta
    }
    fn get_mut_coordinates_patch(&mut self) -> &mut [f64;4]{
        &mut self.coordinates
    }
    fn get_mut_momenta_patch(&mut self) -> &mut [f64;4]{
        &mut self.momenta
    }

    fn get_coordinates_and_momenta(&self) -> [[f64;4];2]{
        [*self.get_coordinates_patch(),*self.get_momenta_patch()]
    }

    fn get_contravariant_momenta (&self) -> [f64;4]{
        self.metric.to_contra( self.get_coordinates_patch(), self.get_momenta_patch() )    
    }

    fn covariant_derivative (& self,index:i8, coor: & [f64;4],  vec :& [f64;4] ) ->  Box<dyn Fn(f64)->f64+ 'a>{
        self.metric.covariant_derivative(index, coor,vec )
    }

    fn get_cartesian_coordinates_and_momenta(&self) -> [f64;8] {
        let coor = self.get_coordinates_patch();
        let mom = self.get_momenta_patch();
        
        [coor[0], coor[1],coor[2],coor[3],
        mom[0], mom[1], mom[2], mom[3]]
    }

    fn sanitize_coordinates(&mut self) {

    }
}


// 

pub fn spher_to_cart(coor: &[f64;3], mom: &[f64;3]) -> [ [f64;3];2] {

    let ps = coor[1].sin();
    let pc = coor[1].cos();
    let hs = coor[2].sin();
    let hc = coor[2].cos();
    
    let pr = mom[0];
    let ph = mom[2];
    let pp = mom[1];
    
    let r = coor[0];
    
    [
        [
            r* hs*pc,
            r* hs*ps,
            r* hc
        ], 
        [
            hs*pc*pr  - r*hs*ps*pp + r*hc*pc*ph,
            hs*ps*pr  + r*hs*pc*pp+ r*hc*ps*ph,
            hc*pr +0.0-r*hs*ph
        ]
    ]

}

pub fn cart_to_spher (coordinates: &[f64;3], momenta: &[f64;3]) -> [ [f64;3];2] {

    let r = ( coordinates[0].powi(2) + coordinates[1].powi(2) + coordinates[2].powi(2)).sqrt();

    let theta = (coordinates[2]/r).acos();
    let phi =  coordinates[1].atan2(coordinates[0]);

    let r2 =  (coordinates[0].powi(2) + coordinates[1].powi(2) ).sqrt();

    // let pc = coordinates[0]/r2;
    // let ps = coordinates[1]/r2;
    // let hc = coordinates[2]/r;
    // let hs = r2 /r;

    let pc = phi.cos();
    let ps = phi.sin();
    let hc = theta.cos();
    let hs = theta.sin();

    let newpos = [
        r,
        phi,
        theta];

    let newmom =  [ 
    hs*pc* momenta[0] + hs*ps*momenta[1]+ hc*momenta[2],
    -ps/(r*hs)*momenta[0] + pc/(r*hs) *momenta[1],
    hc*pc/r *momenta[0] + hc*ps/r*momenta[1]-hs/r *momenta[2]];

    [newpos,newmom]

}