extern crate csv;
#[macro_use]
extern crate serde_derive;

use std::error::Error;
use csv::Writer;
use std::process::Command;

const GLOBAL_METRIC : SchwarzschildMetric = SchwarzschildMetric{r_s:1.0} ;

// metric signature +---
// vector are in x^alpha form   (p_mu p^mu = m^2 c^2)
// for ease of calculation c = hbar = 1
struct SchwarzschildMetric {
   r_s : f64,
}


// christoffel pseudotensor in the form Gamma^(alpha)_(beta,gamma), symmetric in last 2 indices
// coordinates: t,r,theta,phi (= t,r,h,p)
// vanishing symbols are not explicitly defined
#[derive(Debug, Copy, Clone)]
struct SchwarzschildChristoffels {
    trt: f64,
    rrr: f64,  
    rtt: f64, 
    rpp: f64,
    rhh: f64,
    hrh: f64,
    prp: f64,
    hpp: f64,
    php: f64,
}


impl SchwarzschildChristoffels {
    fn update (&mut self, coor : &SchwarzschildCoordinate) {
        let r = coor.r;
        let rs = GLOBAL_METRIC.r_s;
        let a = r-rs;
        let theta = coor.h;
        let st = theta.sin();
        let ct = theta.cos();

        self.trt = rs/ (2.0*r*a);
        self.rrr = -self.trt;
        self.rtt = rs*a/(2.0* r.powi(3));
        self.rpp = -a* st.powi(2);
        self.rhh = -a;
        self.hrh = 1.0/r;
        self.prp = self.hrh;
        self.hpp = -st*ct;
        self.php = ct/st;
     }
     fn new() -> SchwarzschildChristoffels{
        SchwarzschildChristoffels{
            trt: 0.0,
            rrr: 0.0,  
            rtt: 0.0, 
            rpp: 0.0,
            rhh: 0.0,
            hrh: 0.0,
            prp: 0.0,
            hpp: 0.0,
            php: 0.0,
        }
    }
}

// usenumes instead of index for indexing array with state vector, t means x^0, pt dx^0/d tau
#[derive(Copy, Clone,Serialize)]
struct SchwarzschildCoordinate {
    t: f64,
    r: f64,
    h: f64,
    p: f64,
    pt: f64, // = hbar*k
    pr: f64,
    ph: f64,
    pp: f64
} 

impl std::fmt::Debug for SchwarzschildCoordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!( f, "t {:.3},r {:.3},theta {:.3},phi {:.3},pt {:.3},pr {:.3},ptheta {:.3},pphi: {:.3}", self.t ,self.r ,self.h ,self.p ,self.pt, self.pr ,self.ph ,self.pp)
    }
}

#[derive( Copy, Clone)]
struct CartCoordinate {
    t: f64,
    x: f64,
    y: f64,
    z: f64,
    pt: f64, // = hbar*k
    px: f64,
    py: f64,
    pz: f64
} 

impl std::fmt::Debug for CartCoordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!( f, "t {:.3},x {:.3},y {:.3},z {:.3},pt {:.3},px {:.3},py {:.3},pz: {:.3}", self.t ,self.x ,self.y ,self.z ,self.pt, self.px ,self.py ,self.pz)
    }
}

impl std::fmt::Display for CartCoordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!( f, "t {:.3},x {:.3},y {:.3},z {:.3},pt {:.3},px {:.3},py {:.3},pz: {:.3}", self.t ,self.x ,self.y ,self.z ,self.pt, self.px ,self.py ,self.pz)
    }
}

struct Photon {
    state: SchwarzschildCoordinate, //the four generalized coordinates and associated momenta
    christoffels : SchwarzschildChristoffels,
    steps_left: i32,
}

impl std::fmt::Display for Photon {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!( f, "{:?}", self.state)
    }
}

impl Photon {
    fn take_step (&mut self) {
        self.steps_left -=1;

        if self.steps_left == 0 || self.state.r < GLOBAL_METRIC.r_s*1.05 {
            return;
        }

        // update the christoffel symbols based on the current point in space
        self.christoffels.update( &self.state);

        let copy_christoffels = self.christoffels; //creates copy
        let copy_coordinates = self.state;

        //closures for the RK4 integration based on geodesic equation  (Dp^mu/D tau = 0 for force free field)
        let dp_r = |p_r: f64| -> f64 { -(copy_christoffels.rrr* p_r.powi(2) +  copy_christoffels.rtt *  copy_coordinates.pt.powi(2) + copy_christoffels.rpp *  copy_coordinates.pp.powi(2) + copy_christoffels.rhh *  copy_coordinates.ph.powi(2)   )};
        let dp_h = |p_h: f64| -> f64 { -(2.0*copy_christoffels.hrh* copy_coordinates.pr* p_h + copy_christoffels.hpp* copy_coordinates.pp.powi(2))};
        let dp_p = |p_p: f64| -> f64 { -(2.0*copy_christoffels.prp* copy_coordinates.pr* p_p + 2.0*copy_christoffels.php* copy_coordinates.ph* p_p)};
        let dp_t = |p_t: f64| -> f64 { -(2.0*copy_christoffels.trt* copy_coordinates.pr* p_t) };
        
        let d_lambda = 0.01; //lambda is defined as tau/m

        //update the momenta 
        rk4( dp_r, &mut self.state.pr , d_lambda);
        rk4( dp_h, &mut self.state.ph , d_lambda);
        rk4( dp_p, &mut self.state.pp , d_lambda);
        rk4( dp_t, &mut self.state.pt , d_lambda);

        // closures for RK4 integratin   (D x^mu/ D tau = P^mu -> d x_mu/dtau = p^mu + Gamma^mu_alpha,beta x^alpha x^beta
        


        //update coordinates
        rk4( |_| copy_coordinates.pr, &mut self.state.r , d_lambda);
        rk4( |_| copy_coordinates.ph, &mut self.state.h , d_lambda);
        rk4(|_| copy_coordinates.pp, &mut self.state.p , d_lambda);
        rk4( |_| copy_coordinates.pt, &mut self.state.t , d_lambda);

        //todo check for collision with target
    }
}


/// integrate equation dy/dt = f(t,y) over a time step dt and store result in orig params
fn rk4 (f : impl Fn(f64) -> f64, y: &mut f64,dt:f64) {
    let k1 = dt* f(*y);
    let k2 = dt* f(*y+k1/2.0);
    let k3 = dt* f(*y+k2/2.0);
    let k4 = dt* f(*y+k3);
    *y=*y+(k1+2.0*k2+2.0*k3+k4)/6.0;
}


// no correct way to convert t, should be used with caution
fn cart_to_spherical(cart : CartCoordinate ) -> SchwarzschildCoordinate{
    let r = ( cart.x.powi(2) + cart.y.powi(2) + cart.z.powi(2)).sqrt();

    let theta = (cart.z/r).acos();
    let phi =  cart.y.atan2(cart.x);

    let r2 =  (cart.x.powi(2) + cart.y.powi(2) ).sqrt();

    let pc = cart.x/r2;
    let ps = cart.y/r2;
    let hc = cart.z/r;
    let hs = r2 /r;
    
    SchwarzschildCoordinate{
        t: cart.t,
        r: r,
        h: theta,
        p: phi,
        pt: cart.pt, // still E/c = hbar k
        // https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
        pr: hs*pc* cart.px + hs*ps*cart.py+ hc*cart.pz,
        ph: hc*pc/r *cart.px + hc*ps/r*cart.py-hs/r *cart.pz,
        pp: -ps/(r*hs)*cart.px + pc/(r*hs) *cart.py,
    }
}

// 
fn spherical_to_cart(spher : SchwarzschildCoordinate ) -> CartCoordinate{
   let ps = spher.p.sin();
   let pc = spher.p.cos();
   let hs = spher.h.sin();
   let hc = spher.h.cos();

   let pr = spher.pr;
   let ph = spher.ph;
   let pp = spher.pp;

   let r = spher.r;

    CartCoordinate{
        t : spher.t,
        x : r* hs*pc,
        y : r* hs*ps,
        z : r* hc,
        pt: spher.pt,
        // https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
        px: hs*pc*pr + r*hc*pc*ph - r*hs*ps*pp, 
        py: hs*ps*pr + r*hc*ps*ph + r*ps*pc*pp,
        pz: hc*pr -r*hs*ph,
    }
}


// generate photons on line (x,y,z) = (-distance,0,i*spacing ), i=0..number going in the +x direction
fn generate_parallel_photons( distance : f64, spacing: f64, number: i32 ) -> Vec<Photon> {
    let mut a :Vec<Photon> = Vec::new(); 
    
    for i in 0..number {
        a.push( Photon {
            state: cart_to_spherical(CartCoordinate{
                t: 0.0,
                x: -distance,
                y: 0.0,
                z: (i as f64)*spacing,
                pt: 1.0, // = hbar*k
                px: 1.0,
                py: 0.0,
                pz: 0.0,
            }), 
            christoffels : SchwarzschildChristoffels::new(),
            steps_left: 100,
        });
    }
    return a;
}

fn launch_python(n :i32){
    let outp = Command::new("python3")
            .arg("src/plotter.py")
            .arg(format!("-N {}",n) )
            .spawn()
            .expect("error while calling python");

}

fn save_to_csv(v : &Vec<SchwarzschildCoordinate>, name : String) -> Result<(), Box<Error>> {
    let mut wtr = Writer::from_path(  name )?;

    for x in v {
        wtr.serialize( x )?;
    }

    wtr.flush()?;
    Ok(())
}


fn main() {

    // print!( "{}",spherical_to_cart( cart_to_spherical( CartCoordinate{
    //         t: 0.0,
    //         x: -5.0,
    //         y: 3.0,
    //         z: 7.0,
    //         pt: 1.0, // = hbar*k
    //         px: 1.6,
    //         py: 2.90,
    //         pz: 0.35,
    //     })) )  ;


    let numPhotons = 50;
    let mut results : Vec<Vec<SchwarzschildCoordinate>> = Vec::new();


    let mut photons  = generate_parallel_photons(5.0,0.2,numPhotons);
     for (i,x) in photons.iter_mut().enumerate() {

        let mut v : Vec<SchwarzschildCoordinate> = Vec::new();

         // println!("{}",x);
         println!("Begin Pos: {:?}", spherical_to_cart(x.state) );
         for _ in 0..2000 {
            v.push(  x.state  ); //creates a copy to push
            x.take_step();
         } 

         let _ = save_to_csv( &v, format!("files/photon{}.csv", i));

         results.push(v);
     }

     launch_python(numPhotons);
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