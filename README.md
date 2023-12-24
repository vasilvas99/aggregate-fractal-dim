# Calculate fractal dimension of an aggregate

## Usage

```shell
$ aggregate-fractal-dim --help
A CLI tool that takes 3D+t aggregation simulations as 4D *.npz matrices and calculates the fractal dimension of the aggregate

Usage: aggregate-fractal-dim [OPTIONS] <NPZ_FILE_PATH>

Arguments:
  <NPZ_FILE_PATH>  Path to the simulation output

Options:
  -o, --output-file <OUTPUT_FILE>      Path to the output file (CSV) [default: fractal_dimension.csv]
  -s, --csv-separator <CSV_SEPARATOR>  [default: "\t"]
  -h, --help                           Print help
  -V, --version                        Print version
```

## Binaries

Get pre-built binaries from the releases section: [Releases](https://github.com/vasilvas99/aggregate-fractal-dim/releases)
