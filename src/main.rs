pub mod curved_space;
pub mod ray_tracer;
pub mod python_interface;
pub mod test_setups;
pub mod black_body;
pub mod video_maker;




fn main() {
    
    //test_setups::launch_parallel_photons_kerr();

    //look under generated files for output

    //println!("plotting in flat space!");
    //test_setups::ray_trace_minkowski();

    //slower and older code
    //println!("plotting in curved space!");
    //test_setups::ray_trace_schwarzshild();
    
    //println!("plotting in kerr space!");
    //test_setups::ray_trace_kerr();

    video_maker::make_kerr_video();
    //python_interface::make_video(58,"kerr2",10);

    // /old tests

    //test_setups::generate_blackbody(1000.0,10000.0,200);
    //test_setups::photon_orbit();
    //test_setups::test_coordinate_system();
    
}



