extern crate csv;

use std::error::Error;
use std::process::Command;

pub fn launch_python(n :i32, mode : &str){
    Command::new("python3")
            .arg("src/plotter.py")
            .arg(format!("-N {}",n) )
            .arg(format!("-M {}",mode) )
            .spawn()
            .expect("error while calling python");
}

pub fn save_to_csv<S: serde::Serialize>(v : &Vec< S >, name : String) -> Result<(), Box<dyn Error>> {
    let mut wtr = csv::Writer::from_path(  name )?;
    for x in v {
        wtr.serialize( x )?;
    }
    wtr.flush()?;
    Ok(())
}



pub fn make_video(n :i32, base_name: &str, frame_rate:i32){

    Command::new( "ffmpeg")
        .args( &[
        "-framerate",&format!("{}",frame_rate),
        "-i", &format!("video_files/{}%04d.bmp",base_name),
        "-x265-params", "lossless=1", //lossless h265
        "-frames:v", &format!("{}",n),
        "-y", //force overwrite
        &format!("generated_files/output_{}.mkv",base_name),
        ])
        .spawn()
        .expect("failed to execute process");

   
      
}