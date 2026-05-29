# Hybrid Python-C++ Integration for CBN

This repository supports a hybrid execution model where network orchestration and topology generation are handled in Python, while heavy computations (Step 1: Attractors, Step 2: Compatible Pairs, Step 3: Attractor Fields) are offloaded to a high-performance C++ engine.

## Repository Structure

- `cbn_cpp/`: The high-performance C++ engine.
- `cbn_core/`: The Python core submodule (orchestration).
- `python/cbn_engine_utils.py`: Utilities to interface Python with the C++ engine.
- `examples/`: Hybrid execution examples.

## Integration Instructions

To integrate this engine into your thesis repository (`ufabc_phd_cbn`):

1. **Add Submodule**:
   ```bash
   mkdir -p phd_experiments && cd phd_experiments
   git submodule add https://github.com/crptcusco/cbn_complex_networks.git cbn_engine
   git submodule update --init --recursive
   ```

2. **Compile C++ Engine**:
   ```bash
   cd phd_experiments/cbn_engine/cbn_cpp
   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make
   ```

## Hybrid Workflow Example

You can now use `export_cbn_to_json` to bridge Python and C++:

```python
import sys
import os
from cbnetwork.cbnetwork import CBN
from python.cbn_engine_utils import export_cbn_to_json, run_cpp_engine

# 1. Generate network in Python
o_cbn = CBN.cbn_generator(v_topology=1, n_local_networks=4, ...)

# 2. Export and Run C++ Engine
json_path = "network_config.json"
export_cbn_to_json(o_cbn, json_path)
metrics = run_cpp_engine(json_path, binary_path="./cbn_engine/cbn_cpp/build/run_pipeline")

print(f"Attractor Fields: {metrics['metrics']['attractor_fields_count']}")
```

Refer to `examples/hybrid_execution_example.py` and `examples/exp8_data_hybrid.py` for detailed implementations.
