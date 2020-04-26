
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

    //getters

    fn get_mass(&self)->f64;

    fn get_coordinates_patch(&self) -> &[f64;4];
    fn get_momenta_patch(&self) -> &[f64;4];

    fn get_mut_coordinates_patch(&mut self) -> &mut [f64;4];
    fn get_mut_momenta_patch(&mut self) -> &mut [f64;4];

    //calculated properties
    fn calculate_coordinates_and_contravariant_momenta(&self) -> [[f64;4];2]; 
    fn calculate_cartesian_coordinates_and_momenta(&self) -> [[f64;4];2];
    fn calculate_covariant_momenta (&self) -> [f64;4];

    //implementation specific functions
   
    fn sanitize_coordinates(&mut self);
    fn estimate_d_lambda(&self)-> f64;  
    fn get_vector_derivative(&self, coor: &[f64;4], mom: &[f64;4] ) -> [[f64;4];2];
    fn get_error_estimate(&self)->f64;

    // general functions
    fn take_step(&mut self, d_lambda : f64)-> f64 {
        let k0_mom = self.get_momenta_patch();
        let k0_coor = self.get_coordinates_patch();
    
        let rk4_helper= |prev_k : Option< [[f64;4];2]>, factor: f64 | -> [[f64;4];2] {
            let mut knew_arg_coor: [f64;4] = k0_coor.clone();
            let mut knew_arg_mom: [f64;4] = k0_mom.clone();

            match prev_k {
                Some(x)=> {
                    for i in 0..4{
                        knew_arg_mom[i] += factor*x[1][i];
                        knew_arg_coor[i] += factor*x[0][i];
                    }
                },
                None => {},
            }

            let [coor_derivative, mom_derivative] = self.get_vector_derivative( &knew_arg_coor, &knew_arg_mom );
            let mut new_coor = [0.0, 0.0, 0.0, 0.0];
            let mut new_mom = [0.0, 0.0, 0.0, 0.0];

            for i in 0..4{
                new_coor[i] = coor_derivative[i]*d_lambda;
                new_mom[i] = mom_derivative[i]*d_lambda;
            }

            [new_coor, new_mom]
        };
 
        // here rk4 scheme
        let [k1_coor,k1_mom] = rk4_helper( None, 0.0 );
        let [k2_coor,k2_mom] = rk4_helper( Some( [k1_coor,k1_mom] ), 0.5 );
        let [k3_coor,k3_mom] = rk4_helper( Some( [k2_coor,k2_mom] ), 0.5 );
        let [k4_coor,k4_mom] = rk4_helper( Some( [k3_coor,k3_mom] ), 1.0 );

        // synthesis
        let mut_coor = self.get_mut_coordinates_patch();
        for i in 0..4{
            mut_coor[i] =  mut_coor[i] + (k1_coor[i] + 2.0*k2_coor[i] + 2.0*k3_coor[i]+k4_coor[i] )/6.0;
        }

        let mut_mom = self.get_mut_momenta_patch();
        for i in 0..4{
            mut_mom[i] = mut_mom[i] + (k1_mom[i] + 2.0*k2_mom[i] + 2.0*k3_mom[i]+k4_mom[i] )/6.0;
        }

        return self.get_error_estimate();
    }

    fn restore_coordinates(&mut self, coordinates: [f64;4], momenta: [f64;4], patch: u8){
        let coor = self.get_mut_coordinates_patch();
        *coor = coordinates;

        let mom = self.get_mut_momenta_patch();
        *mom = momenta;

        self.set_patch(patch);
        self.sanitize_coordinates();
    }

    fn get_patch(&self)->u8;
    fn set_patch(&mut self, patch: u8);

    // debug
    fn print(&self){
        let [coor,mom] = self.calculate_coordinates_and_contravariant_momenta();
        println!("coord: [{:.5},{:.5},{:.5},{:.5}] mom:[{:.5},{:.5},{:.5},{:.5}]", coor[0],coor[1],coor[2],coor[3],mom[0],mom[1],mom[2],mom[3])
    }

}

pub trait Metric<'a> : std::marker::Sync + std::marker::Send {
    type SpecificSpaceObject : SpaceObject<'a> + std::marker::Sync + std::marker::Send ;

    //fn get_christoffel (&self,x0: u8,x1: u8,x2: u8,coor: &[f64;4]) -> f64;
    fn g_lower(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64;
    fn g_upper(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64;

    fn to_covariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4];
    fn to_contravariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4];

    fn spawn_space_object(&'a self, coordinates : [f64;4], momenta : [f64;4] , mass : f64 ) -> Self::SpecificSpaceObject;
    fn spawn_space_object_from_cartesian(&'a self, coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject;

}

////////////////////implementation of the schwarzschildmetric

//t,r,phi,theta
pub struct SchwarzschildMetric {
    pub r_s : f64,
    pub delta: f64,
    pub max_step : f64,
}

pub fn new_schwarzschild_metric(r_s:f64, delta:f64, max_step: f64) -> SchwarzschildMetric {
    SchwarzschildMetric{r_s: r_s, delta:delta, max_step: max_step}
}

impl<'a> Metric<'a> for SchwarzschildMetric {
    type SpecificSpaceObject= SchwarzschildObject<'a>;

    

    fn g_lower(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64{
        match (x0,x1) {
            (0,0) => (1.0 - self.r_s/coor[ 1 ] ),
            (1,1) => -1.0/(1.0 - self.r_s/coor[ 1 ] ),
            (2,2) => -(coor[1] * coor[3].sin()  ).powi(2),
            (3,3) => -coor[1].powi(2),
            __ => 0.0,
        }
    }

    fn g_upper(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64{
        1.0/self.g_lower(x0,x1,coor) 
    }

    fn to_covariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ self.g_lower( 0, 0, coordinates )* vec[0], self.g_lower( 1, 1, coordinates )* vec[1]  ,self.g_lower( 2, 2, coordinates )* vec[2],self.g_lower( 3, 3, coordinates )* vec[3] ]
    }

    fn to_contravariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ self.g_upper( 0, 0, coordinates )* vec[0], self.g_upper( 1, 1, coordinates )* vec[1]  ,self.g_upper( 2, 2, coordinates )* vec[2],self.g_upper( 3, 3, coordinates )* vec[3] ]
    }

    fn spawn_space_object(&'a self ,coordinates : [f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{ 
        let mut  res = SchwarzschildObject{  
             patch_coordinates : coordinates,
             patch_momenta: momenta,
             mass: mass, metric: &self, 
             coordinate_patch: 0, 
             //theta_threshold: 0.5*std::f64::consts::PI/2.0, 
             killing_const: [ 0.0, 0.0 ],
        };        

        res.sanitize_coordinates();
        res.reset_killing_const();

        res
    }

    fn spawn_space_object_from_cartesian( &'a self ,coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{
        let [coor,mom] = cart_to_spher(   &[coordinates[1],coordinates[2],coordinates[3]] ,  &[momenta[1],momenta[2],momenta[3]] );

        let newpos = [
            coordinates[0],
            coor[0],
            coor[1],
            coor[2]
        ];
        let p_t = ( mass.powi(2) -(mom[0].powi(2)*self.g_lower(1,1,&newpos) +mom[1].powi(2)*self.g_lower(2,2,&newpos)+mom[2].powi(2)*self.g_lower(3,3,&newpos))/( self.g_lower(0,0,&newpos) )  ).sqrt();

        let newmom =  [ 
            p_t, 
            mom[0],
            mom[1],
            mom[2]
        ];

        self.spawn_space_object( newpos , newmom, mass)
    }

    
}

impl SchwarzschildMetric {

    // D x^mu/(D tau) = d x^mu/(d tau) + gamma^mu_alpha,beta x^alpha x^beta, this returns a closure |x^mu| gamma^mu_alpha,beta x^alpha x^beta 
    fn covariant_derivative (&self, coordinates: &[f64;4], vector :& [f64;4]  ) ->   [f64;4]{
        let vect = vector.clone();
        let coor = coordinates.clone();        
        [
            {
                let trt = self.get_christoffel(0,1,0, &coor);
                -(2.0* trt *  vect[1]* vect[0])
            },
            {
                let rrr = self.get_christoffel(1,1,1, &coor);
                let rtt = self.get_christoffel(1,0,0, &coor);
                let rpp = self.get_christoffel(1,2,2, &coor);
                let rhh = self.get_christoffel(1,3,3, &coor);

                -(rrr* vect[1].powi(2) +  rtt *  vect[0].powi(2) + rpp *  vect[2].powi(2) + rhh *  vect[3].powi(2)   )
            },
            {
                let prp = self.get_christoffel(2,1,2, &coor);
                let php = self.get_christoffel(2,3,2, &coor);

                -(2.0*prp* vect[1]* vect[2] + 2.0*php* vect[3]* vect[2])
            },
            {
                let hrh = self.get_christoffel(3,1,3, &coor);
                let hpp = self.get_christoffel(3,2,2, &coor);
                -(2.0*hrh* vect[1]* vect[3] + hpp* vect[2].powi(2))
            },
        ]
    } 

    fn get_christoffel (&self,x0: u8,x1: u8,x2: u8,coor: &[f64;4]) -> f64{
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

}

//implementation of a object living in a schwarzschildmetric

pub struct SchwarzschildObject <'a> {
    patch_coordinates : [f64;4],
    patch_momenta :  [f64;4],
    metric: &'a SchwarzschildMetric,
    mass : f64,
    coordinate_patch: u8,
    killing_const: [f64;2],
}

impl<'a> SpaceObject<'a> for SchwarzschildObject<'a>{

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

    fn calculate_coordinates_and_contravariant_momenta(&self) -> [[f64;4];2] {
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

    fn calculate_covariant_momenta (&self) -> [f64;4]{
        self.metric.to_covariant_vector( self.get_coordinates_patch(), self.get_momenta_patch() )    
    }

    fn get_vector_derivative(&self, coor: &[f64;4], mom: &[f64;4] ) -> [[f64;4];2]{
        return [ mom.clone(), self.metric.covariant_derivative( coor, mom)  ];
    }

    fn calculate_cartesian_coordinates_and_momenta(&self) -> [[f64;4];2] {
        let mut ret : [[f64;4];2] = [[0.0;4];2];
        let [coordinates,momenta] = self.calculate_coordinates_and_contravariant_momenta();

        let [coor,mom] = spher_to_cart(  &[coordinates[1],coordinates[2],coordinates[3]] ,  &[momenta[1],momenta[2],momenta[3]]);

        

        ret[0][0] = coordinates[0];
        ret[1][0] = momenta[0];

        for i in 1..4 {
            ret[0][i] = coor[i-1];
            ret[1][i] = mom[i-1];
        }

        ret
    }

    fn sanitize_coordinates(&mut self) {
        
        let mut mut_coor = self.patch_coordinates;
        let mut mut_mom = self.patch_momenta;

        if mut_coor[3] < 0.0 {
            mut_coor[3]   *= -1.0 ;
            mut_mom[3] *= -1.0;
            mut_coor[2] -= std::f64::consts::PI;

        } else if mut_coor[3] >std::f64::consts::PI {
            mut_coor[3]  = 2.0*std::f64::consts::PI -  mut_coor[3];
            mut_mom[3] *= -1.0;
            mut_coor[2] += std::f64::consts::PI;

            //self.reset_killing_const();
        }

        //self.patch_coordinates[2] = self.patch_coordinates[2]%(std::f64::consts::PI*2.0);

        //x,y,z -> yzx
        // if  (std::f64::consts::PI/2.0 - self.patch_coordinates[3]).abs() > self.theta_threshold {
            
        //     //self.switch_patch();
        // }
    }

    fn get_mass(&self)->f64{
        self.mass
    }



    fn get_patch(&self)->u8{
        self.coordinate_patch
    }

    fn set_patch(&mut self, patch: u8){
        self.coordinate_patch = patch;
    }

    //got estimate from https://arxiv.org/pdf/1303.5057.pdf

    fn get_error_estimate(&self)->f64{
        let metr = self.metric;
        let mom = self.patch_momenta;
        let coor =  self.get_coordinates_patch();

        let xi=(metr.g_lower(1,1,coor)* mom[1].powi(2)+
                metr.g_lower(2,2, coor)* mom[2].powi(2)+
                metr.g_lower(3,3, coor)* mom[3].powi(2))/
               (metr.g_lower(0,0, coor)* mom[0].powi(2));

        return (xi +1.0).abs()

    }

    //got estimate from https://arxiv.org/pdf/1303.5057.pdf

    fn estimate_d_lambda(&self)-> f64{
        let mom = self.patch_momenta;
        let coor = self.patch_coordinates;

        let estimate1 = self.metric.delta/( (mom[1]/coor[1]).abs()+mom[3].abs()+mom[2].abs());
        let estimate2 = (coor[1]-self.metric.r_s)/(2.0* mom[1].abs());
        
        let min = if estimate1 < estimate2 {
            estimate1
        }else{

            //println!("r estimate");
            estimate2
        };

        return if self.metric.max_step< min {
           self.metric.max_step
        } else {
            min
        }  
    }
}

//XYZ->YZX

impl<'a> SchwarzschildObject<'a>{

    // fn update_momenta_killing_symmetry(&mut self) -> f64 {
    //     //in schwarzschild metric: not explicitly dependend upon: t and phi =>d/d lambda (k_t^mu p_mu) =0 and   d/d lambda (k_phi^mu p_mu) =0 
    //     let coor = self.get_coordinates_patch().clone();

    //     let g00 = self.metric.g_upper(0,0, &coor );
    //     let g22 = self.metric.g_upper(2,2, &coor);

    //     let [e,l] = self.killing_const.clone(); // E,L,C

       
    //     let mut_mom = self.get_mut_momenta_patch();

    //     let c0 = e*g00;
    //     let c2 = l*g22;

       
    //     mut_mom[0] = c0;
    //     mut_mom[2] = c2;


    //     0.0
    // }


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

        self.reset_killing_const();
    }

    fn reset_killing_const(&mut self){
        let mom = self.get_momenta_patch();
        
        let e = mom[0]*self.metric.g_lower(0,0,self.get_coordinates_patch());
        let l = mom[2]*self.metric.g_lower(2,2,self.get_coordinates_patch());

        self.killing_const = [ e ,l];
    }
}


////////////////////implementation of minkowski metric 




//t,x,y,z
pub struct MinkowskiMetric {
    d_lambda:f64,
}

pub fn new_minkowski_metric(d_lambda:f64) -> MinkowskiMetric {
    MinkowskiMetric{ d_lambda: d_lambda }
}


impl<'a> Metric<'a> for MinkowskiMetric {
    type SpecificSpaceObject= MinkowskiObject<'a>;

    fn g_lower(&self, x0 :u8, x1 :u8,_ :&[f64;4]) -> f64{
        match (x0,x1) {
            (0,0) => 1.0 ,
            (1,1) => -1.0,
            (2,2) => -1.0,
            (3,3) => -1.0,
            __ => 0.0,
        }
    }

    fn g_upper(&self, x0 :u8, x1 :u8,_ :&[f64;4]) -> f64{
        match (x0,x1) {
            (0,0) => 1.0 ,
            (1,1) => -1.0,
            (2,2) => -1.0,
            (3,3) => -1.0,
            __ => 0.0,
        }
    }

    fn to_covariant_vector(&self,_ : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ vec[0], - vec[1]  ,- vec[2],- vec[3]]
    }

    fn to_contravariant_vector(&self,_ : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ vec[0], - vec[1]  ,- vec[2],- vec[3]]
    }


    fn spawn_space_object(&'a self ,coordinates : [f64;4], momenta : [f64;4] , mass : f64 ) -> Self::SpecificSpaceObject{ 
        MinkowskiObject{ coordinates : coordinates, momenta: momenta, mass: mass, metric: &self }
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
    mass : f64,
}

impl<'a> SpaceObject<'a> for MinkowskiObject<'a>{

    fn get_mass(&self)->f64{
        self.mass
    }
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

    fn calculate_coordinates_and_contravariant_momenta(&self) -> [[f64;4];2]{
        [*self.get_coordinates_patch(),*self.get_momenta_patch()]
    }

    fn calculate_covariant_momenta (&self) -> [f64;4]{
        self.metric.to_covariant_vector( self.get_coordinates_patch(), self.get_momenta_patch() )    
    }

    fn get_vector_derivative(&self, _: &[f64;4], mom: &[f64;4] ) -> [[f64;4];2] {
        return [ mom.clone(), [0.0,0.0,0.0,0.0]  ];
    }

    fn calculate_cartesian_coordinates_and_momenta(&self) -> [[f64;4];2] {
        let coor = self.get_coordinates_patch();
        let mom = self.get_momenta_patch(); 
        [ coor.clone(),mom.clone() ]
    }

    fn sanitize_coordinates(&mut self) {
    }

    fn get_patch(&self)->u8{
       0
    }

    fn set_patch(&mut self, _: u8){
    }

    fn get_error_estimate(&self)->f64{
        0.0
    }

    fn estimate_d_lambda(&self)->f64{
        self.metric.d_lambda
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

    let pc = coordinates[0]/r2;
    let ps = coordinates[1]/r2;
    let hc = coordinates[2]/r;
    let hs = r2 /r;


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