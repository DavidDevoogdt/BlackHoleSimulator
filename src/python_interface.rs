extern crate csv;

use std::error::Error;
use std::process::Command;

pub fn launch_python(n :i32){
    Command::new("python3")
            .arg("src/plotter.py")
            .arg(format!("-N {}",n) )
            .spawn()
            .expect("error while calling python");

}

pub fn save_to_csv(v : &Vec< [f64;8] >, name : String) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::Writer::from_path(  name )?;

    for x in v {
        wtr.serialize( x )?;
    }

    wtr.flush()?;
    Ok(())
}
