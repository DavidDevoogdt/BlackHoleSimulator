[![](http://img.youtube.com/vi/a6b87AsqNdc/0.jpg)](http://www.youtube.com/watch?v=a6b87AsqNdc "Black Hole vido")


# Kerr Black Hole Simulator
## installation
install the rust toolchain. [info here](https://www.rust-lang.org/tools/install)

to run it, go to the root folder of this project, and run:
```
cargo run --release
```
In src/main.rs, the different components can be activated,  
in src/test_setups.rs, the parameters of the simulation can be tweaked,
## overview
This is hobby project which implements a general relativistic ray tracer in rust.
  
The core of the raytracer is in src/curved_space.rs. Here, a trait (class) is provided for different metrics togheter with a trait for the objects living in the corresponding space. The numerical stepping procedures, conversion between contra and covariant vectors, spawning objects from cartesian representation etc are implemented here.  

src/ray_tracer takes as input objects such as an accretion disk and a background picture and a discription of the camera object. The light rays are propageted until they collide and then the color is deduced. It also provides methods to generate output image and plots in python. 

src/test_setup is a predefined set of input argumets for the ray_tracer.

## internals
The rays are propagated with a runge kutta stepping sheme.
  
For the schwarzschild metric, the contraviriant momenta are updated with the gedesic equition (with christffel symbols). It uses 2 coordinate patches to avoid difficulties at the poles. The stepping is done with RK4 and an specific error estimate and a adaptive step size based on the current coordinates and momenta.
  
For the kerr metric, the covariant momenta are updated using embedded rk5 algoritm and the update equation from [this](https://arxiv.org/abs/1601.02063) paper.


## pictures
In kerr spacetime, the simultion needs about 40s for a full HD picture on a i7-4700MQ CPU @ 2.40GHz (released in 2013).

Some example pictures:
![Image description](/pictures/kerr-40s.bmp)
![Image description](/pictures/kerr.bmp)
![Image description](/pictures/example_plot.png)
