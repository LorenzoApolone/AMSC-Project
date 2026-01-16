import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

def plot_xy(x, y, x_name="x", y_name="y", save_path="plot.png"):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title(f"{y_name} vs {x_name}")
    plt.grid(True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_result_file(file_path, x_col, y_col, output_plot_path):
    with open(file_path, "r") as f:
        header_line = f.readline().strip()

    column_names = header_line.split()

    if x_col >= len(column_names) or y_col >= len(column_names):
        raise ValueError("x_col or y_col is out of range.")

    x_name = column_names[x_col]
    y_name = column_names[y_col]

    data = np.loadtxt(file_path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    x = data[:, x_col]
    y = data[:, y_col]

    plot_xy(x, y, x_name=x_name, y_name=y_name, save_path=output_plot_path)


# -----------------------------------------------------------------------------
# Updated Columns and Directory Parameters
# -----------------------------------------------------------------------------

# File now includes the following columns:
#   ["max_iter", "tol", "it_n", "delta_x", "final_t"]
FILE_COLUMNS = ["max_iter", "tol", "it_n", "delta_x", "final_t"]

# Directory structure NO LONGER includes stop_idx
# It is now ONLY:
#   test_f / dim / n_points / n_cores /
DIR_PARAMS = ["test_f", "dim", "n_points", "n_cores"]


def load_file_data(file_path):
    with open(file_path, "r") as f:
        header = f.readline().strip().split()
    data = np.loadtxt(file_path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return header, data


def extract_dir_params(path_parts):
    """Extract directory parameters from path parts (NO stop_idx)."""
    return {
        "test_f": path_parts[-4],
        "dim": int(path_parts[-3]),
        "n_points": int(path_parts[-2]),
        "n_cores": int(path_parts[-1]),
    }


def satisfies_fixed_params(dir_params, file_data, fixed_params):
    header, data = file_data

    for key, val in fixed_params.items():

        # -----------------------------
        # Parameters coming from folders
        # -----------------------------
        if key in DIR_PARAMS:
            if dir_params[key] != val:
                return False

        # -----------------------------
        # Parameters coming from file columns
        # -----------------------------
        elif key in FILE_COLUMNS:
            idx = header.index(key)
            if not np.allclose(data[:, idx], val):
                return False

        else:
            raise ValueError(f"Unknown fixed parameter: {key}")

    return True


def extract_x(dir_params, file_data, x):
    header, data = file_data

    if x in DIR_PARAMS:
        return np.array([float(dir_params[x])] * len(data))

    if x in FILE_COLUMNS:
        return data[:, header.index(x)]

    raise ValueError(f"Unknown x parameter: {x}")


def extract_y(dir_params, file_data, y):
    header, data = file_data

    if y in DIR_PARAMS:
        return np.array([float(dir_params[y])] * len(data))

    if y in FILE_COLUMNS:
        return data[:, header.index(y)]

    raise ValueError(f"Unknown y parameter: {y}")


# -----------------------------------------------------------------------------
# Aggregated Plots
# -----------------------------------------------------------------------------

def plot_average_results(fixed_params, x, y, save_path):
    tests_dir = Path("./tests")
    curves = []

    for test_f_dir in tests_dir.iterdir():
        if not test_f_dir.is_dir(): continue

        for dim_dir in test_f_dir.iterdir():
            if not dim_dir.is_dir(): continue

            for np_dir in dim_dir.iterdir():
                if not np_dir.is_dir(): continue

                for nc_dir in np_dir.iterdir():
                    if not nc_dir.is_dir(): continue

                    dir_params = extract_dir_params(nc_dir.parts)

                    for file in nc_dir.glob("*.txt"):
                        header, data = load_file_data(file)
                        file_data = (header, data)

                        if not satisfies_fixed_params(dir_params, file_data, fixed_params):
                            continue

                        x_arr = extract_x(dir_params, file_data, x)
                        y_arr = extract_y(dir_params, file_data, y)

                        curves.append((x_arr, y_arr))

    # If nothing matched, avoid empty plot
    if not curves:
        print(f"[WARNING] No matching runs found for {fixed_params}. Empty plot skipped.")
        return

    bucket = defaultdict(list)
    for x_arr, y_arr in curves:
        for xi, yi in zip(x_arr, y_arr):
            bucket[xi].append(yi)

    xs = np.array(sorted(bucket.keys()))
    means = np.array([np.mean(bucket[xv]) for xv in xs])
    stds  = np.array([np.std(bucket[xv]) for xv in xs])

    upper = means + 2 * stds
    lower = means - 2 * stds

    plt.figure(figsize=(10, 6))
    plt.plot(xs, means, label="mean y(x)")
    plt.fill_between(xs, lower, upper, alpha=0.3, label="±2σ band")

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{y} vs {x} (mean ± 2σ)")
    plt.grid(True)
    plt.legend()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_median_results(fixed_params, x, y, save_path):
    tests_dir = Path("./tests")
    curves = []

    for test_f_dir in tests_dir.iterdir():
        if not test_f_dir.is_dir(): continue

        for dim_dir in test_f_dir.iterdir():
            if not dim_dir.is_dir(): continue

            for np_dir in dim_dir.iterdir():
                if not np_dir.is_dir(): continue

                for nc_dir in np_dir.iterdir():
                    if not nc_dir.is_dir(): continue

                    dir_params = extract_dir_params(nc_dir.parts)

                    for file in nc_dir.glob("*.txt"):
                        header, data = load_file_data(file)
                        file_data = (header, data)

                        if not satisfies_fixed_params(dir_params, file_data, fixed_params):
                            continue

                        x_arr = extract_x(dir_params, file_data, x)
                        y_arr = extract_y(dir_params, file_data, y)

                        curves.append((x_arr, y_arr))

    if not curves:
        print(f"[WARNING] No matching runs found for {fixed_params}. Empty plot skipped.")
        return

    bucket = defaultdict(list)
    for x_arr, y_arr in curves:
        for xi, yi in zip(x_arr, y_arr):
            bucket[xi].append(yi)

    xs = np.array(sorted(bucket.keys()))
    medians = np.array([np.median(bucket[xv]) for xv in xs])

    plt.figure(figsize=(10, 6))
    plt.plot(xs, medians, label="median y(x)")

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(f"{y} vs {x} (median)")
    plt.grid(True)
    plt.legend()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate plots from swarm search experiment outputs."
    )

    # -------------------------------------------------------------------------
    # Choice of plot type
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--mode",
        type=str,
        default="average",
        choices=["average", "median"],
        help="Type of aggregate plot to generate.",
    )

    # -------------------------------------------------------------------------
    # x and y axis
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--x",
        type=str,
        required=True,
        help=f"x-axis parameter (one of {DIR_PARAMS + FILE_COLUMNS})",
    )

    parser.add_argument(
        "--y",
        type=str,
        required=True,
        help=f"y-axis parameter (one of {DIR_PARAMS + FILE_COLUMNS})",
    )

    # -------------------------------------------------------------------------
    # Fixed parameters
    # Allow: test_f, dim, n_points, n_cores, max_iter, tol, it_n, delta_x, final_t
    # -------------------------------------------------------------------------
    for p in DIR_PARAMS + FILE_COLUMNS:
        parser.add_argument(
            f"--{p}",
            type=float,
            help=f"Fix the parameter '{p}' to a specific value.",
        )

    # -------------------------------------------------------------------------
    # Output file
    # -------------------------------------------------------------------------
    parser.add_argument(
        "--out",
        type=str,
        required=True,
        help="Output file path for the generated plot (.png recommended).",
    )

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # Build fixed_params dictionary
    # Only include parameters the user actually provided
    # -------------------------------------------------------------------------
    fixed_params = {}
    for p in DIR_PARAMS + FILE_COLUMNS:
        val = getattr(args, p)
        if val is not None:
            # Convert integer-valued directory params properly
            if p in ["dim", "n_points", "n_cores"]:
                val = int(val)
            fixed_params[p] = val

    # -------------------------------------------------------------------------
    # Dispatch plot type
    # -------------------------------------------------------------------------
    if args.mode == "average":
        plot_average_results(fixed_params, args.x, args.y, args.out)
    else:
        plot_median_results(fixed_params, args.x, args.y, args.out)
