import sys
import os
import time
import pandas as pd

# Add project root and cbn_core/src to path
sys.path.append(os.path.abspath("."))
sys.path.append(os.path.abspath("cbn_core/src"))

from cbnetwork.cbnetwork import CBN
from cbnetwork.globaltopology import GlobalTopology
from python.cbn_engine_utils import export_cbn_to_json, run_cpp_engine

def run_exp8_hybrid():
    # Simulation Parameters
    n_nets_list = [3, 4]
    n_vars = 6
    topologies = [1, 3] # Complete, Cycle

    results = []

    print(f"{'Topology':<15} | {'N_Nets':<6} | {'Status':<10} | {'Total Time (s)':<15}")
    print("-" * 55)

    for top_id in topologies:
        top_name = GlobalTopology.allowed_topologies[top_id]
        for n_nets in n_nets_list:
            # 1. Generate System
            o_cbn = CBN.cbn_generator(
                v_topology=top_id,
                n_local_networks=n_nets,
                n_vars_network=n_vars,
                n_input_variables=2,
                n_output_variables=2,
                n_max_of_clauses=3,
                n_max_of_literals=2
            )

            # 2. Hybrid Execution
            json_path = f"tmp_net_{top_id}_{n_nets}.json"
            export_cbn_to_json(o_cbn, json_path)

            try:
                metrics = run_cpp_engine(json_path)
                m = metrics["metrics"]
                print(f"{top_name:<15} | {n_nets:<6} | {metrics['status']:<10} | {m['total_time']:>15.4f}")

                results.append({
                    "topology": top_name,
                    "n_networks": n_nets,
                    "n_vars": n_vars,
                    "step1_time": m["step1_time"],
                    "step2_time": m["step2_time"],
                    "step3_time": m["step3_time"],
                    "total_time": m["total_time"],
                    "attractors": m["local_attractors_count"],
                    "fields": m["attractor_fields_count"]
                })
            except Exception as e:
                print(f"Error in {top_name} ({n_nets}): {e}")
            finally:
                if os.path.exists(json_path):
                    os.remove(json_path)

    # 3. Export Consolidated Results
    df = pd.DataFrame(results)
    df.to_csv("exp8_hybrid_results.csv", index=False)
    print("\nResults saved to exp8_hybrid_results.csv")

if __name__ == "__main__":
    run_exp8_hybrid()
