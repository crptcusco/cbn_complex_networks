# CBN C++ High-Performance Engine

This directory contains the pure C++ implementation of the CBN simulation pipeline, optimized for mass data production in HPC environments.

## Features
- **Native Generation**: Stochastic topology and dynamics generation (Complete, Cycle, Path digraphs).
- **HPC Optimized**: Constant memory footprint via strict scoping and pre-allocated vectors.
- **Resilient Logging**: Direct CSV output with immediate flushing.
- **SLURM Ready**: CLI support for easy integration with Job Arrays.

## Compilation

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Running Experiments

The `experiment_runner` binary is the main entry point for data production.

```bash
./experiment_runner --samples 1000 --networks 5 --vars 10 --topology 1 --output_file results_complete_5_10.csv
```

### Options:
- `--samples <int>`: Number of systems to generate and analyze.
- `--topology <int>`: 1 (Complete), 3 (Cycle), 4 (Path).
- `--networks <int>`: Number of local networks.
- `--vars <int>`: Variables per local network.
- `--inputs <int>`: Input signals per local network.
- `--outputs <int>`: Output variables per local network.
- `--output_file <string>`: Destination CSV.

## Diagnostic Tool

Use `run_pipeline` to execute the algorithm on a specific JSON configuration:

```bash
./run_pipeline network_config.json
```
