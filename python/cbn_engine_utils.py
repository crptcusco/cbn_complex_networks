import json
import subprocess
import os

def export_cbn_to_json(o_cbn, filepath):
    """
    Serializes a CBN object to a JSON file compatible with the C++ engine.
    """
    data = {
        "local_networks": [],
        "directed_edges": []
    }

    for net in o_cbn.l_local_networks:
        net_data = {
            "index": net.index,
            "internal_variables": net.internal_variables,
            "descriptive_function_variables": []
        }
        for var in net.descriptive_function_variables:
            net_data["descriptive_function_variables"].append({
                "index": var.index,
                "cnf": var.cnf_function
            })
        data["local_networks"].append(net_data)

    for edge in o_cbn.l_directed_edges:
        edge_data = {
            "index": edge.index,
            "index_variable": edge.index_variable,
            "input_local_network": edge.input_local_network,
            "output_local_network": edge.output_local_network,
            "output_variables": edge.l_output_variables,
            "coupling_function": edge.coupling_function,
            "true_table": edge.true_table
        }
        data["directed_edges"].append(edge_data)

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=4)

def run_cpp_engine(json_path, binary_path="./cbn_cpp/build/run_pipeline"):
    """
    Invokes the C++ engine and returns the parsed metrics.
    """
    if not os.path.exists(binary_path):
        raise FileNotFoundError(f"C++ binary not found at {binary_path}")

    result = subprocess.run([binary_path, json_path], capture_output=True, text=True, check=True)
    return json.loads(result.stdout)
