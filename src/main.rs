use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};
use clap::Parser;
use fractal_analysis::*;
use ndarray::{Array4, ArrayView3};
use npyz::Deserialize;

static ARR_DEFAULT_NAME: &str = "arr_0";
static MAX_VAL: i32 = 255;

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

    let mut npz = npyz::npz::NpzArchive::open(file_path)?;
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

fn threshold(x: &i32) -> i32 {
    if *x < 2 {
        return 0;
    }
    return MAX_VAL;
}

fn calculate_fractal_dimension_3d(frame: ArrayView3<i32>) -> f64 {
    let frame = frame.map(threshold);
    let s = frame.shape();
    let w = s[0];
    let l = s[1];
    let h = s[2];

    let get_key_from_sample = |(coords, val): ((usize, usize, usize), &i32)| -> u32 {
        let (x, y, z) = coords;
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, w, 0);
        let norm_y = |y| normalise_as_u8(y, l, 0);
        let norm_z = |z| normalise_as_u8(z, h, 0);

        // Safety: Ensure values are clipped to 255;
        let arr = [norm_x(x), norm_y(y), norm_z(z), *val as u8];
        morton_encoding::morton_encode(arr)
    };

    let clzs = get_clzs(frame.indexed_iter(), get_key_from_sample);
    let (tmp, lacun) = get_results_from_clzs(clzs);
    let res = finalise_results::<32>(tmp, lacun, w * l * h, 8);
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
        .quote_style(csv::QuoteStyle::NonNumeric)
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
