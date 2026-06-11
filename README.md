# Coupled Boolean Network (CBN) Engine - C++ Implementation

This repository contains a high-performance C++ engine for simulating and benchmarking Coupled Boolean Networks (CBNs). It provides robust implementations of various execution strategies, ranging from traditional sequential processing to advanced parallel algorithms using OpenMP and asynchronous task reduction.

## Core Features

- **High-Performance Logic:** Optimized truth table generation and Boolean update function evaluation.
- **Advanced Parallelism:**
  - `TraditionalExperiment`: Single-threaded baseline for performance comparison.
  - `SimpleParallelExperiment`: OpenMP-accelerated phase execution.
  - `AdvancedParallelExperiment`: Features asynchronous binary reduction and bucket-balanced workload distribution for maximum scalability.
- **Flexible Topologies:** Supports various network topologies (Complete, Cycle, Scale-Free, etc.) via `cbn_generator`.
- **Scientific Benchmarking:** Integrated pipeline for large-scale performance data collection.

## Compilation Instructions

The project uses CMake and requires a C++17 compliant compiler (e.g., GCC 7+, Clang 5+) and OpenMP.

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

This will generate several executables, including the primary benchmarking tool: `scientific_benchmarking`.

## Scientific Benchmarking Tool

The `scientific_benchmarking` utility is designed for rigorous performance evaluation across multiple randomized samples and strategies.

### Usage

```bash
./scientific_benchmarking [OPTIONS]
```

### CLI Arguments

- `--samples <int>`: Number of independent randomized samples to evaluate (default: 3).
- `--networks` or `--nodes <int>`: Number of local networks in the CBN (default: 6).
- `--vars <int>`: Number of internal Boolean variables per local network (default: 10).
- `--topology <int>`: Global topology type (default: 1).
- `--output <path>`: Filename for the global metrics CSV (default: `benchmark_results.csv`).
- `--dir <path>`: Directory where detailed sample data (fields and topology) will be stored (default: `phd_experiments/output_data`).

### Example

To run 100 samples of a 11-node CBN with topology 4:

```bash
./scientific_benchmarking --samples 100 --nodes 11 --vars 10 --topology 4 --output results_topo4.csv
```

## Output Data Structure

The benchmarking tool generates three types of persistence artifacts:

1.  **Global Log (`benchmark_results.csv`):** A cumulative CSV file containing timing metrics (Phase 1, 2, and 3), memory usage (RSS), and attractor counts for every sample and strategy execution.
2.  **Attractor Fields JSON (`cbn_sample_{ID}_{STRATEGY}_fields.json`):** Contains the set of stable attractor fields found for a specific sample and strategy.
3.  **Topology CSV (`cbn_sample_{ID}_{STRATEGY}_topology.csv`):** Exports the directed edges of the generated CBN.

All detailed files are saved in the directory specified by `--dir`.

## Architecture Overview

- `cbnetwork/`: Core logic for Boolean networks, edges, and topologies.
- `experiment_strategies.hpp`: Strategy pattern implementation for different execution paradigms.
- `examples/`: Utility programs and benchmarking drivers.

## License

This project is part of doctoral research at UFABC. All rights reserved.
