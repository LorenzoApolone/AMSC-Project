import json
import numpy as np
from pathlib import Path
from plots import plot_average_results, plot_median_results
from collections import defaultdict

# ------------------------------
# Load parameter combinations
# ------------------------------
param_file = Path("./run_parameters.json")
if not param_file.exists():
    raise FileNotFoundError(f"{param_file} not found. Run the C++ script first.")

with open(param_file) as f:
    param_list = json.load(f)

plots_base = Path("./plots")
plots_base.mkdir(exist_ok=True)


# ------------------------------
# Mapping agg_type -> (default x, default y, fixed param keys)
# ------------------------------
agg_mapping = {
    0: {"x": "it_n", "y": "delta_x", "fixed_keys": ["dim", "n_points", "max_iter", "tol"]},
    1: {"x": "n_cores", "y": "final_t", "fixed_keys": ["dim", "n_points", "max_iter", "tol"]},
    2: {"x": "dim", "y": "final_t", "fixed_keys": ["n_points", "n_cores", "max_iter", "tol"]},
    3: {"x": "dim", "y": "it_n", "fixed_keys": ["n_points", "max_iter", "tol"]},
}

# ------------------------------
# Scan tests directory to see which parameter combinations exist
# ------------------------------
tests_dir = Path("./tests")
existing_runs = []

for test_f_dir in tests_dir.iterdir():
    if not test_f_dir.is_dir():
        continue

    # Extract function name
    function_name = test_f_dir.name

    # Make plots/<function_name>
    plots_dir = plots_base / function_name
    plots_dir.mkdir(exist_ok=True)

    for dim_dir in test_f_dir.iterdir():
        if not dim_dir.is_dir():
            continue

        for npoints_dir in dim_dir.iterdir():
            if not npoints_dir.is_dir():
                continue

            for ncores_dir in npoints_dir.iterdir():
                if not ncores_dir.is_dir():
                    continue

                # (rest of existing logic)


                # Find all files inside this folder
                for file in ncores_dir.glob("*.txt"):

                    # Read header + first data row
                    with open(file, "r") as f:
                        header = f.readline().strip().split()
                    data = np.loadtxt(file, skiprows=1)
                    if data.ndim == 1:
                        data = data.reshape(1, -1)

                    # Extract max_iter and tol
                    max_iter = float(data[0, header.index("max_iter")]) if "max_iter" in header else None
                    tol      = float(data[0, header.index("tol")])      if "tol" in header else None

                    # Append full configuration
                    existing_runs.append({
                        "test_f": test_f_dir.name,
                        "dim": int(dim_dir.name),
                        "n_points": int(npoints_dir.name),
                        "n_cores": int(ncores_dir.name),
                        "max_iter": max_iter,
                        "tol": tol,
                    })


if not existing_runs:
    raise RuntimeError("No runs found in ./tests!")

# ------------------------------
# Loop over existing runs and agg_types
# ------------------------------
for run_params in existing_runs:
    for agg_type, mapping in agg_mapping.items():
        # Build fixed_params dictionary from keys in mapping
        fixed_params = {k: run_params[k] for k in mapping["fixed_keys"] if k in run_params}

        # Construct save file name
        save_file_name = f"agg_type_{agg_type}_" + "_".join(f"{k}{v}" for k, v in fixed_params.items()) + ".png"
        save_file = plots_dir / save_file_name

        # Skip if plot already exists
        if save_file.exists():
            continue

        print(f"Generating plot: {save_file}")
        x_val = mapping["x"]
        y_val = mapping["y"]

        # Call the plotting function (mean aggregation)
        plot_average_results(fixed_params, x_val, y_val, save_file)

print("All existing-run plots generated successfully.")
