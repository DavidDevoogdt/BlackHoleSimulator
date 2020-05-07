extern crate csv;
extern crate image;
extern crate serde;
use self::serde::{Deserialize, Serialize};

use std::fmt;

pub struct CieLookup {
    start: i32,
    end: i32,
    delta: i32,
    steps: usize,
    array: Vec<XYZ>, // array of [X,Y,Z] structs
}

#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug)]
struct XYZ {
    lambda: i32,
    X: f64,
    Y: f64,
    Z: f64,
}

impl std::fmt::Display for XYZ {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}nm, [{:.04},{:.04},{:.04}]",
            self.lambda, self.X, self.Y, self.Z
        )
    }
}

impl CieLookup {
    pub fn new() -> CieLookup {
        let mut arr: Vec<XYZ> = vec![];

        // Build the CSV reader and iterate over each record.
        //data table taken from: http://cvrl.ucl.ac.uk/cmfs.htm
        let mut rdr = csv::Reader::from_path("src_files/CIE_color_table.csv").unwrap();

        for result in rdr.records() {
            let record = result.unwrap();
            let row: XYZ = record.deserialize(None).unwrap();

            arr.push(row);
        }

        let ret = CieLookup {
            start: arr[0].lambda,
            steps: arr.len(),
            end: arr.last().expect("something wen wrong parsing").lambda,
            delta: arr[1].lambda - arr[0].lambda,
            array: arr,
        };

        return ret;
    }

    pub fn get_black_body_intensity(&self, temp: f64, wavelength_shift: f64) -> image::Rgb<u8> {
        // https://en.wikipedia.org/wiki/Planckian_locus
        //T in kelvin, lambda in nm
        let planck = |lambda: f64, temp: f64| -> f64 {
            5.0 * 1.0E15 * lambda.powi(-5) / ((1.43877E7 / (temp * lambda)).exp() - 1.0)
        };

        let mut res_x: f64 = 0.0;
        let mut res_y: f64 = 0.0;
        let mut res_z: f64 = 0.0;

        for i in 0..self.steps {
            let xyz = &self.array[i];

            let pl = planck(xyz.lambda as f64 + wavelength_shift, temp);
            res_x += pl * xyz.X;
            res_y += pl * xyz.Y;
            res_z += pl * xyz.Z;
        }

        //normalize for case of finer grained integrtion
        res_x *= self.delta as f64;
        res_y *= self.delta as f64;
        res_z *= self.delta as f64;

        Self::xyz_to_rgb(res_x, res_y, res_z)
    }

    pub fn wavelength_to_rgb(&self, lambda: f64) -> image::Rgb<u8> {
        if lambda < self.start as f64 || lambda > self.end as f64 {
            return image::Rgb([0, 0, 0]);
        }

        let low = (lambda / (self.delta as f64)).floor() as i32 - self.start / self.delta;

        let xyz1 = &self.array[low as usize];
        let xyz2 = &self.array[(low + 1) as usize];

        let frac: f64 = (lambda - xyz1.lambda as f64) / (self.delta as f64);

        return Self::xyz_to_rgb(
            frac * xyz1.X + (1.0 - frac) * xyz2.X,
            frac * xyz1.Y + (1.0 - frac) * xyz2.Y,
            frac * xyz1.Z + (1.0 - frac) * xyz2.Z,
        );
    }

    pub fn xyz_to_rgb(x: f64, y: f64, z: f64) -> image::Rgb<u8> {
        //https://en.wikipedia.org/wiki/CIE_1931_color_space
        let mut big_r = 0.41847 * x + -0.15866 * y + -0.082835 * z;
        let mut big_g = -0.091169 * x + 0.25243 * y + 0.015708 * z;
        let mut big_b = 0.00092090 * x + -0.0025498 * y + 0.17860 * z;

        if big_r < 0.0 {
            big_r = 0.0;
        }
        if big_g < 0.0 {
            big_g = 0.0;
        }
        if big_b < 0.0 {
            big_b = 0.0;
        }

        println!("{} {} {}", big_r, big_g, big_b);

        let r: u8 = (255.0 * big_r / (big_r + big_g + big_b)) as u8;
        let g: u8 = (255.0 * big_g / (big_r + big_g + big_b)) as u8;
        let b: u8 = (255.0 * big_b / (big_r + big_g + big_b)) as u8;

        println!("{} {} {}", r, g, b);

        return image::Rgb([r, g, b]);
    }
}
