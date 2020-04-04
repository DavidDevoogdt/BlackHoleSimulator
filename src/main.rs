pub mod curved_space;
pub mod ray_tracer;
pub mod python_interface;
pub mod test_setups;


fn main() {
    test_setups::launch_parallel_photons();
    test_setups::ray_trace();
 }

