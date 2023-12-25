use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};
use clap::Parser;
use fractal_analysis::*;
use ndarray::{Array4, ArrayView3};
use npyz::Deserialize;
use rayon::prelude::*;

static ARR_DEFAULT_NAME: &str = "arr_0";

/// A CLI tool that takes 3D+t aggregation simulations
/// as 4D *.npz matrices and calculates the fractal dimension
/// of the aggregate.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Path to the simulation output
    #[arg()]
    npz_file_path: PathBuf,

    /// Path to the output file (CSV)
    #[arg(short = 'o', long, default_value = "fractal_dimension.csv")]
    output_file: PathBuf,

    #[arg(short = 's', long, default_value_t = '\t')]
    csv_separator: char,
}

#[derive(Debug, serde::Serialize)]
#[serde(rename_all = "PascalCase")]
struct CsvRecord {
    frame_number: usize,
    fractal_dimension: f64,
}

fn load_aggregate_data<T: Deserialize>(file_path: impl AsRef<Path>) -> Result<Array4<T>> {
    use ndarray::ShapeBuilder;
    let file = std::io::Cursor::new(std::fs::read(file_path)?); // Read the whole file in one shot
    let mut npz = npyz::npz::NpzArchive::new(file)?;
    let arr = npz
        .by_name(ARR_DEFAULT_NAME)?
        .ok_or_else(|| anyhow!("Could not load array by name {}", ARR_DEFAULT_NAME))?;
    let shape = arr.shape().to_vec();
    let order = arr.order();
    let data: Vec<T> = arr.into_vec()?;

    let shape = match shape[..] {
        [i1, i2, i3, i4] => [i1 as usize, i2 as usize, i3 as usize, i4 as usize],
        _ => return Err(anyhow!("expected 4D array")),
    };
    let true_shape = shape.set_f(order == npyz::Order::Fortran);

    Ok(ndarray::Array4::from_shape_vec(true_shape, data)?)
}

fn threshold(x: &i32) -> u8 {
    if *x < 2 {
        return u8::MIN;
    }
    return u8::MAX;
}

fn calculate_fractal_dimension_3d(frame: ArrayView3<i32>) -> f64 {
    let frame = frame.map(threshold);
    let s = frame.shape();
    let x_max = s[0];
    let y_max = s[1];
    let z_max = s[2];
    let buf = frame.into_raw_vec().into_par_iter().enumerate();

    let get_key_from_sample = |(flattened_coord, val): (usize, u8)| -> u32 {
        let mut idx = flattened_coord;
        let z = idx / (x_max * y_max);
        idx -= z * x_max * y_max;
        let y = idx / x_max;
        let x = idx % x_max;

        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, x_max, 0);
        let norm_y = |y| normalise_as_u8(y, y_max, 0);
        let norm_z = |z| normalise_as_u8(z, z_max, 0);

        let arr = [norm_x(x), norm_y(y), norm_z(z), val];
        lindel::morton_encode(arr)
    };

    let clzs = get_clzs_par(buf, get_key_from_sample).collect::<Vec<_>>();
    let (tmp, lacun) = get_results_from_clzs(clzs.into_iter());
    let res = finalise_results::<32>(tmp, lacun, x_max * y_max * z_max, 8);
    res.0
}

fn get_separator(sep_char: char) -> Result<u8> {
    let mut buf = [0; 1];
    sep_char.encode_utf8(&mut buf);

    Ok(buf[0])
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let p = PathBuf::from(cli.npz_file_path);
    let l = load_aggregate_data::<i32>(p)?;
    println!("Loading done. Starting processing.");

    let output_file = std::fs::File::create(cli.output_file)?;
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(get_separator(cli.csv_separator)?)
        .from_writer(output_file);

    for (frame_number, frame) in l.outer_iter().enumerate() {
        let fractal_dimension = calculate_fractal_dimension_3d(frame);
        wtr.serialize(CsvRecord {
            frame_number,
            fractal_dimension,
        })?;
        if frame_number % 10 == 0 {
            wtr.flush()?;
        }
        println!("Processed frame: {frame_number}");
    }

    wtr.flush()?;

    Ok(())
}
