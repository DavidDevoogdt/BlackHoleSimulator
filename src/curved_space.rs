
////////////////////////////////////////////////////////////////////
// trait definitions
////////////////////////////////////////////////////////////////////


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
    
    fn get_d_lambda(&self)-> f64;
    fn set_d_lambda(&mut self, d_lambda: Option<f64>);
    
    fn get_vector_derivative(&self, coor: &[f64;4], mom: &[f64;4] ) -> [[f64;4];2];

    fn get_error_estimate(&self)->f64;
    fn set_error_estimate(&mut self, err: Option<f64>);

    // Implement stepping mechanism or chose from default implementations below
    fn take_step(&mut self, d_lambda : f64);

    // general step functions
    // rk4 stepping with external error estimate
    fn rk4_stepper(&mut self, d_lambda : f64) {
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

        self.set_error_estimate(None);
        self.set_d_lambda(None);
        
    }

    // rk4 scheme with error estimation and stepsize  (embedded  method)
    // for numbers, see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
    fn rk5_stepper(&mut self, d_lambda : f64, precision: f64) {
        let butcher_tableau: [[f64;5];5] = [
            [1.0/4.0,          0.0,                0.0,            0.0,            0.0],
            [3.0/32.0,         9.0/32.0,           0.0,            0.0,            0.0],
            [1932.0/2197.0,    -7200.0/2197.0,     7296.0/2197.0,  0.0,            0.0],
            [439.0/216.0,      -8.0,     	        3680.0/513.0,   -845.0/4104.0,  0.0],
            [-8.0/27.0, 	    2.0,      	        -3544.0/2565.0, 1859.0/4104.0,  -11.0/40.0],
        ];

        let _ = [ 
            0.0,
            1.0/4.0,
            3.0/8.0,
            12.0/13.0,
            1.0,
            0.5
        ];

        let b_coeff = [
            16.0/135.0,	
            0.0, 
            6656.0/12825.0,	
            28561.0/56430.0, 
            -9.0/50.0,	
            2.0/55.0
        ];

        let b_star_coeff = [
            25.0/216.0,	
            0.0,	
            1408.0/2565.0,	
            2197.0/4104.0,	
            -1.0/5.0,	
            0.0,
        ];


        let y0_mom = self.get_momenta_patch();
        let y0_coor = self.get_coordinates_patch();

        let mut k_i: Vec< [[f64;4];2] > = Vec::with_capacity(6);

        for i in 0..6 {
            let [coor_arg,mom_arg ] = {
                let [mut new_coor,mut new_mom] = [
                    [0.0,0.0,0.0,0.0],
                    [0.0,0.0,0.0,0.0],
                ];

                for j in 0..i{
                    let fact = butcher_tableau[i-1][j];
                    for s in 0..4 {
                        new_coor[s] += fact* k_i[j][0][s];
                        new_mom[s] += fact* k_i[j][1][s];
                    }
                }

                for i in 1..4{
                    new_coor[i] += y0_coor[i];
                    new_mom[i] += y0_mom[i];
                }

                [new_coor, new_mom]
            };

            let mut new_k_i = self.get_vector_derivative( &coor_arg, &mom_arg  );

            for i in 0..4{
                new_k_i[0][i] *= d_lambda;
                new_k_i[1][i] *= d_lambda;
            }

            k_i.push( new_k_i );
        }

        

        // synthesis
        let mut_coor = self.get_mut_coordinates_patch();
       
        for j in 0..6{
            for i in 0..4{
                mut_coor[i] += k_i[j][0][i]*b_star_coeff[j];
            }
        }

        let mut_mom = self.get_mut_momenta_patch();
        for j in 0..6{
            for i in 0..4{
                mut_mom[i] +=  k_i[j][1][i]*b_star_coeff[j];
            }
        }

        let err: f64 = {
            let [mut err_coor,mut err_mom] = [
                [0.0,0.0,0.0,0.0],
                [0.0,0.0,0.0,0.0],
            ];

            for j in 0..6{
                let fact = b_coeff[j]-b_star_coeff[j];
                for i in 1..4 {
                    err_coor[i] += fact*k_i[j][0][i];
                    err_mom[i] += fact*k_i[j][1][i];
                }

            }

            let mut err = 0.0;
            for i in 0..4{
                err += err_coor[i].powi(2) + err_mom[i].powi(2);
            }

            ( err/8.0 ).sqrt()
        };

        self.set_error_estimate(Some(err));
        self.set_d_lambda( Some(0.7*d_lambda*  ( precision/err ).powf(0.2) )  );//0.9: safety factor

    }

    // used to recover from inaccurate steps
    fn restore_state(&mut self);
    fn store_state(&mut self);

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

////////////////////////////////////////////////////////////////////
// implemntation of schwarzschild metric
////////////////////////////////////////////////////////////////////

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
             error_estimate: 0.0,

             backup_patch_coordinates : coordinates,
             backup_patch_momenta: momenta, 
             backup_coordinate_patch: 0, 
             d_lambda: 0.0,
        };        

        res.sanitize_coordinates();
        res.reset_killing_const();

        res.store_state();

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



pub struct SchwarzschildObject <'a> {
    patch_coordinates : [f64;4],
    patch_momenta :  [f64;4],
    coordinate_patch: u8,
    metric: &'a SchwarzschildMetric,
    mass : f64,
    killing_const: [f64;2],
    error_estimate: f64,

    backup_patch_coordinates : [f64;4],
    backup_patch_momenta :  [f64;4],
    backup_coordinate_patch: u8,
    d_lambda: f64,
}

impl<'a> SpaceObject<'a> for SchwarzschildObject<'a>{

    fn take_step(&mut self, d_lambda : f64){
        self.rk4_stepper(d_lambda);
    }

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

        self.patch_coordinates[2] = self.patch_coordinates[2]%(std::f64::consts::PI*2.0);

        //x,y,z -> yzx
        if  (std::f64::consts::PI/2.0 - self.patch_coordinates[3]).abs() > std::f64::consts::PI/4.0 {
            
            self.switch_patch();
        }
    }

    fn get_mass(&self)->f64{
        self.mass
    }

    fn restore_state(&mut self){
        self.patch_coordinates = self.backup_patch_coordinates.clone();
        self.patch_momenta= self.backup_patch_momenta.clone();
        self.coordinate_patch = self.backup_coordinate_patch.clone();
    }


    fn store_state(&mut self){
        self.backup_patch_coordinates = self.patch_coordinates.clone();
        self.backup_patch_momenta= self.patch_momenta.clone();
        self.backup_coordinate_patch = self.coordinate_patch.clone();
    }


    fn get_error_estimate(&self)->f64{
        self.error_estimate
    }

    //got estimate from https://arxiv.org/pdf/1303.5057.pdf
    fn set_error_estimate(&mut self,error: Option<f64>){
        match error{
            Some(x) => self.error_estimate =x,
            None => {
                let metr = self.metric;
                let mom = self.patch_momenta;
                let coor =  self.get_coordinates_patch();
        
                let xi=(metr.g_lower(1,1,coor)* mom[1].powi(2)+
                        metr.g_lower(2,2, coor)* mom[2].powi(2)+
                        metr.g_lower(3,3, coor)* mom[3].powi(2))/
                       (metr.g_lower(0,0, coor)* mom[0].powi(2));
        
                self.error_estimate =  (xi +1.0).abs();
            }
        }
       
    }


    fn get_d_lambda(&self) -> f64 {
        return self.d_lambda;
    }

    //got estimate from https://arxiv.org/pdf/1303.5057.pdf
    fn set_d_lambda(&mut self, d_lambda: Option<f64>){
        match d_lambda {
            Some(x) => self.d_lambda = x,
            None=>{
                let mom = self.patch_momenta;
                let coor = self.patch_coordinates;
        
                let estimate1 = self.metric.delta/( (mom[1]/coor[1]).abs()+mom[3].abs()+mom[2].abs());
                let estimate2 = (coor[1]-self.metric.r_s)/(2.0* mom[1].abs());
                
                let min = if estimate1 < estimate2 {
                    estimate1
                }else{
                    estimate2
                };
        
                return if self.metric.max_step< min {
                   self.d_lambda =  self.metric.max_step;
                } else {
                    self.d_lambda = min;
                }  
            }
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


////////////////////////////////////////////////////////////////////
// implemntation of minkowski metric
////////////////////////////////////////////////////////////////////

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
        MinkowskiObject{ coordinates : coordinates, momenta: momenta, mass: mass, metric: &self, error_estimate:0.0,backup_coordinates:[ 0.0, 0.0, 0.0, 0.0], backup_momenta:[ 0.0, 0.0, 0.0, 0.0] }
    }

    fn spawn_space_object_from_cartesian( &'a self ,coordinates :[f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{
        self.spawn_space_object( coordinates , momenta, mass)
    }

   
}

pub struct MinkowskiObject <'a> {
    coordinates : [f64;4],
    momenta :  [f64;4],
    metric: &'a MinkowskiMetric,
    mass : f64,
    error_estimate: f64,
    backup_coordinates : [f64;4],
    backup_momenta :  [f64;4],
}

impl<'a> SpaceObject<'a> for MinkowskiObject<'a>{

    fn take_step(&mut self, d_lambda : f64){
        self.rk4_stepper(d_lambda);
    }

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

   

    fn get_error_estimate(&self)->f64{
        self.error_estimate
    }

    fn set_error_estimate(&mut self, err: Option<f64>){
        match err{
            Some(x) => {
                self.error_estimate = x
            },
            None => panic!("no automatic error estimate for minkowski implemented")
        }
    }

    
    fn get_d_lambda(&self)->f64{
        self.metric.d_lambda
    }

    fn set_d_lambda(&mut self, _: Option<f64>) {
    }

    fn restore_state(&mut self){
        self.coordinates = self.backup_coordinates.clone();
        self.momenta= self.backup_momenta.clone();
      
    }

    fn store_state(&mut self){
        self.backup_coordinates = self.coordinates.clone();
        self.backup_momenta= self.momenta.clone();
    }
}


////////////////////////////////////////////////////////////////////
// implemntation of kerr metric with hamiltonian equations of motion
////////////////////////////////////////////////////////////////////

#[allow(non_snake_case)]
pub struct KerrMetric {
    pub M : f64,
    pub J : f64,
    pub a:f64,
    pub r_s:f64,
    pub step_precision:f64,
}

pub fn new_kerr_metric( j:f64, step_parameter: f64) -> KerrMetric {
    KerrMetric {
        M:1.0,
        J:j,
        a: j,
        r_s:2.0,
        step_precision: step_parameter,
    }
}

impl KerrMetric {

    fn get_sigma(&self, coor :&[f64;4])->f64{
        coor[1].powi(2) + (self.a* coor[2].cos()).powi(2) 
    }

    fn get_delta(&self, coor:&[f64;4] )->f64{
        coor[1].powi(2) - self.r_s * coor[1]+ self.a.powi(2)
    }
}

impl<'a> Metric<'a> for KerrMetric {
    type SpecificSpaceObject= KerrObject<'a>;
    fn g_lower(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64{
        match (x0,x1){
            (0,0)=> -(1.0- self.r_s*coor[1]/self.get_sigma( coor )),
            (1,1)=> self.get_sigma(coor)/self.get_delta(coor),
            (2,2)=> ( coor[1].powi(2)+ self.a.powi(2) *(1.0+ self.r_s*coor[1]/self.get_sigma( coor ) * coor[3].sin().powi(2))   )*coor[3].sin().powi(2),
            (3,3)=> self.get_sigma(coor),
            (0,2)=>  -2.0*self.r_s*coor[1]*self.a* coor[3].sin().powi(2)/self.get_sigma(coor),
            _ => panic!("wrong call to g_lower kerr"),
        }
    }

    fn g_upper(&self, x0 :u8, x1 :u8,coor :&[f64;4]) -> f64{
       
        let rrs = self.r_s*coor[1]/self.get_sigma( coor );    
        let det = coor[3].sin().powi(2)*( self.a.powi(2)*( coor[3].cos().powi(2)*rrs-1.0) + coor[1].powi(2)*(rrs-1.0)    );

        match (x0,x1){
            (0,0)=> self.g_lower(2,2,coor)/det ,
            (1,1)=> self.get_delta(coor)/self.get_sigma(coor) ,
            (2,2)=> self.g_lower(0,0,coor)/det,
            (3,3)=> 1.0/(self.get_sigma(coor) ),
            (0,2)=> -self.g_lower(0,2,coor)/det,
            _ => panic!("wrong call to g_upper kerr"),
        }
    }

    fn to_covariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ 
            vec[0]* self.g_lower(0,0, coordinates) + 0.5*vec[2]* self.g_lower(0,2, coordinates), 
            vec[1]* self.g_lower(1,1, coordinates), 
            vec[2]* self.g_lower(2,2, coordinates) + 0.5*vec[0]* self.g_lower(0,2, coordinates),
            vec[3]* self.g_lower(3,3, coordinates)  
        ]
    }

    fn to_contravariant_vector(&self, coordinates : &[f64;4] ,vec: &[f64;4]) -> [f64;4]{
        [ 
            vec[0]* self.g_upper(0,0, coordinates) +0.5*vec[2]* self.g_upper(0,2, coordinates),
            vec[1]* self.g_upper(1,1, coordinates),
            vec[2]* self.g_upper(2,2, coordinates) +0.5*vec[0]* self.g_upper(0,2, coordinates),
            vec[3]* self.g_upper(3,3, coordinates)  
        ]
    }

    //momenta input are contravariant
    fn spawn_space_object(&'a self ,coordinates : [f64;4], momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{ 
        let mut  res = KerrObject{  
             coordinates : coordinates,
             momenta: self.to_covariant_vector( &coordinates, &momenta),
             mass: mass, metric: &self, 
             constants_of_motion: [ 0.0, 0.0, 0.0 ],
             error_estimate : 0.0,
             backup_coordinates : [ 0.0, 0.0, 0.0, 0.0],
             backup_momenta: [ 0.0, 0.0, 0.0, 0.0],
             d_lambda: 0.01,
             backup_d_lambda: 0.01,
        };        

        res.sanitize_coordinates();
        res.reset_killing_const();

        res
    }

    fn spawn_space_object_from_cartesian( &'a self ,coordinates :[f64;4], contravariant_momenta : [f64;4] , mass : f64) -> Self::SpecificSpaceObject{
        let [coor,mom] = cart_to_kerr(   &[coordinates[1],coordinates[2],coordinates[3]] ,  &[contravariant_momenta[1],contravariant_momenta[2],contravariant_momenta[3]], self.a );

        let newpos = [
            coordinates[0],
            coor[0],
            coor[1],
            coor[2]
        ];

        let a = self.g_lower(0,0,&newpos);
        let b = self.g_lower(0,2,&newpos)* mom[1];
        let c = (mom[0].powi(2)*self.g_lower(1,1,&newpos) + mom[1].powi(2)*self.g_lower(2,2,&newpos)+mom[2].powi(2)*self.g_lower(3,3,&newpos)) - mass.powi(2);

        let p_t = (-b + ( b.powi(2)-4.0*a*c ).sqrt())/(2.0*a);

        let newmom =  [ 
            p_t, 
            mom[0],
            mom[1],
            mom[2]
        ];

        self.spawn_space_object( newpos , newmom, mass)

    }

    
}


pub struct KerrObject <'a> {
    coordinates : [f64;4], // x^mu
    momenta :  [f64;4], //p_\mu = g_\mu \nu p^\nu
    metric: &'a KerrMetric,
    mass : f64,
    constants_of_motion: [f64;3],
    error_estimate : f64,
    backup_coordinates : [f64;4], 
    backup_momenta :  [f64;4],
    d_lambda: f64,
    backup_d_lambda : f64,
}

impl<'a> SpaceObject<'a> for KerrObject<'a>{

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

    fn calculate_coordinates_and_contravariant_momenta(&self) -> [[f64;4];2] {
        [ 
            self.coordinates,
            self.metric.to_contravariant_vector( &self.coordinates, &self.momenta ) 
        ]
    }

    fn calculate_covariant_momenta (&self) -> [f64;4]{
        self.momenta    
    }

    // information taken from https://arxiv.org/pdf/1601.02063.pdf
    #[allow(non_snake_case)]
    fn get_vector_derivative(&self, coor: &[f64;4], mom: &[f64;4] ) -> [[f64;4];2]{
        let mom_contra = self.metric.to_contravariant_vector( coor, mom);
        let metr = self.metric;

        let [E,L,Kappa] = self.constants_of_motion;

        let sin_theta  = coor[3].sin();
        let cos_theta = coor[3].cos();
        
        [
          [
            mom_contra[0] ,
            mom_contra[1],
            mom_contra[2],
            mom_contra[3],
          ],
          [
            0.0,
            1.0/( metr.get_sigma(coor)*metr.get_delta(coor) )*(
                -Kappa*( coor[1]-1.0 )
                + 2.0*coor[1]*( coor[1].powi(2)+ metr.a.powi(2))*E.powi(2)
                -2.0*metr.a*E*L
            ) -2.0*mom[1].powi(2)*( coor[1]-1.0)/metr.get_sigma(coor),
            0.0,
            sin_theta*cos_theta/metr.get_sigma(coor)*(  (L/sin_theta.powi(2)).powi(2) - (self.metric.a*E).powi(2)   ),
          ]
        ]
    }

    fn calculate_cartesian_coordinates_and_momenta(&self) -> [[f64;4];2] {
        let mut ret : [[f64;4];2] = [[0.0;4];2];

        let [coordinates,momenta] = self.calculate_coordinates_and_contravariant_momenta();

        let [coor,mom] = kerr_to_cart(  &[coordinates[1],coordinates[2],coordinates[3]] ,  &[momenta[1],momenta[2],momenta[3]], self.metric.a );

        ret[0][0] = coordinates[0];
        ret[1][0] = momenta[0];

        for i in 1..4 {
            ret[0][i] = coor[i-1];
            ret[1][i] = mom[i-1];
        }

        ret
    }

    fn sanitize_coordinates(&mut self) {
        
        let mut mut_coor = self.coordinates;
        let mut mut_mom = self.momenta;

        if mut_coor[3] < 0.0 {
             mut_coor[3]   *= -1.0 ;
             mut_mom[3] *= -1.0;
             mut_coor[2] -= std::f64::consts::PI;

        } else if mut_coor[3] >std::f64::consts::PI {
            mut_coor[3]  = 2.0*std::f64::consts::PI -  mut_coor[3];
            mut_mom[3] *= -1.0;
            mut_coor[2] += std::f64::consts::PI;
        }

    }

    fn get_mass(&self)->f64{
        self.mass
    }


    fn get_error_estimate(&self)->f64{
        self.error_estimate
    }

    fn set_error_estimate(&mut self, err: Option<f64>){
        match err{
            Some(x) => self.error_estimate = x,
            None => panic!("no automatic error estimate for kerr implemented")
        }
    }


    fn set_d_lambda(&mut self, d_lambda: Option<f64>) {
        match d_lambda {
            Some(x) => self.d_lambda = x,
            None => {
                panic!{"not implemented"};
            },
        }
    }

    fn get_d_lambda(&self)-> f64{
       self.d_lambda
    }

    fn take_step(&mut self, d_lambda : f64){
        self.rk5_stepper(d_lambda, self.metric.step_precision);
    }

    fn restore_state(&mut self){

        self.coordinates = self.backup_coordinates.clone();
        self.momenta= self.backup_momenta.clone();
        self.d_lambda = self.backup_d_lambda;
      
    }

    fn store_state(&mut self){
        self.backup_coordinates = self.coordinates.clone();
        self.backup_momenta= self.momenta.clone();
        self.backup_d_lambda = self.d_lambda;
    }
}

impl<'a> KerrObject<'a>{

    //https://arxiv.org/pdf/1601.02063.pdf

    #[allow(non_snake_case)]
    fn reset_killing_const(&mut self){
        let mom = self.momenta;
        let coor = self.coordinates;
        
        let E = -mom[0];
        let L = mom[2];
        let kappa = mom[3].powi(2) + (L/ coor[3].sin()).powi(2) + (self.metric.a * E ).powi(2); 
        
        self.constants_of_motion = [ E ,L, kappa];
    }
}



pub fn kerr_to_cart(coor: &[f64;3], mom: &[f64;3],a:f64) -> [ [f64;3];2] {

    let ps = coor[1].sin();
    let pc = coor[1].cos();
    let hs = coor[2].sin();
    let hc = coor[2].cos();
    
    let pr = mom[0];
    let ph = mom[2];
    let pp = mom[1];
    
    let r = coor[0];
    
    let a_r = (r.powi(2)+a.powi(2)).sqrt();

    [
        [
            a_r* hs*pc,
            a_r* hs*ps,
            r* hc
        ], 
        [
            r/a_r*hs*pc*pr  - a_r*hs*ps*pp + a_r*hc*pc*ph,
            r/a_r*hs*ps*pr  + a_r*hs*pc*pp+ a_r*hc*ps*ph,
            hc*pr +0.0-r*hs*ph
        ]
    ]


}

pub fn cart_to_kerr (coordinates: &[f64;3], momenta: &[f64;3],a:f64) -> [ [f64;3];2] {

    let phi = coordinates[1].atan2( coordinates[0]);
    let b =  coordinates[0].powi(2) + coordinates[1].powi(2) + coordinates[2].powi(2) - a.powi(2);
    let r = ( (  b + ( b.powi(2) + 4.0* (a*coordinates[2]).powi(2) ).sqrt()   )/2.0 ).sqrt();
    let theta =  (coordinates[2]/r).acos();

    let newpos = [
        r,
        phi,
        theta
    ];

    let a_r = (r.powi(2)+a.powi(2)).sqrt();
    let n = (theta.cos()*a).powi(2)+r.powi(2);

    let hc = theta.cos();
    let hs = theta.sin();
    let pc = phi.cos();
    let ps = phi.sin();
    
    let newmom =  [ 
        a_r*r/n*hs*pc* momenta[0] + a_r*r/n*hs*ps*momenta[1]+ a_r.powi(2)/n*hc*momenta[2],
        -ps/(a_r* hs)*momenta[0] + pc/(a_r* hs) *momenta[1],
        a_r/n*hc*pc*momenta[0] + a_r/n*hc*ps*momenta[1]-r/n*hs*momenta[2]
    ];

    [newpos,newmom]
}
