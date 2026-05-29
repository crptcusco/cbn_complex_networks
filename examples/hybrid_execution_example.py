import sys
import os
import json

# Add project root and cbn_core/src to path
sys.path.append(os.path.abspath("."))
sys.path.append(os.path.abspath("cbn_core/src"))

from cbnetwork.cbnetwork import CBN
from python.cbn_engine_utils import export_cbn_to_json, run_cpp_engine

def main():
    # 1. Generate Synthetic Topology in Python
    print("Generating network in Python...")
    o_cbn = CBN.cbn_generator(
        v_topology=1, # Complete Graph
        n_local_networks=3,
        n_vars_network=5,
        n_input_variables=1,
        n_output_variables=1,
        n_max_of_clauses=2,
        n_max_of_literals=2
    )

    # 2. Export to JSON
    json_path = "temp_network.json"
    export_cbn_to_json(o_cbn, json_path)
    print(f"Network exported to {json_path}")

    # 3. Invoke C++ Engine
    cpp_binary = "./cbn_cpp/build/run_pipeline"
    print("Executing C++ engine...")
    try:
        metrics = run_cpp_engine(json_path, cpp_binary)

        # 4. Integrate results
        print("\n--- Hybrid Execution Metrics ---")
        print(json.dumps(metrics, indent=4))

        if metrics["status"] == "success":
            print("\nIntegration successful!")

    except Exception as e:
        print(f"ERROR: {e}")
    finally:
        if os.path.exists(json_path):
            os.remove(json_path)

if __name__ == "__main__":
    main()
