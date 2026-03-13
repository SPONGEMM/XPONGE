"""
Convert Xponge JSON analysis outputs to CSV files.
"""
import csv
import json
import os


def _ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _write_rows(path, header, rows):
    _ensure_parent(path)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def _write_xy_series(data, output):
    x = data.get("x", [])
    y = data.get("y", [])
    if len(x) != len(y):
        raise ValueError("JSON x/y arrays must have the same length.")
    _write_rows(output, ["x", "y"], zip(x, y))
    stats = data.get("stats")
    if isinstance(stats, dict) and stats:
        stats_output = output[:-4] + "_stats.csv" if output.endswith(".csv") else output + "_stats.csv"
        _write_rows(stats_output, ["key", "value"], [(key, value) for key, value in stats.items()])
    hbonds = data.get("hbonds")
    if isinstance(hbonds, list) and hbonds:
        width = max(len(row) if isinstance(row, list) else 1 for row in hbonds)
        normalized = []
        for row in hbonds:
            if isinstance(row, list):
                normalized.append(row + [""] * (width - len(row)))
            else:
                normalized.append([row] + [""] * (width - 1))
        hbonds_output = output[:-4] + "_hbonds.csv" if output.endswith(".csv") else output + "_hbonds.csv"
        _write_rows(hbonds_output, [f"col_{i}" for i in range(width)], normalized)


def _write_times_files(data, output):
    times = data.get("times", [])
    files = data.get("files", [])
    if len(times) != len(files):
        raise ValueError("JSON times/files arrays must have the same length.")
    _write_rows(output, ["time_ns", "file"], zip(times, files))


def _write_fes(data, output):
    x_centers = data.get("x_centers", [])
    y_centers = data.get("y_centers", [])
    free_energy = data.get("free_energy", [])
    if len(free_energy) != len(x_centers):
        raise ValueError("free_energy row count must match x_centers length.")
    rows = []
    for i, x_val in enumerate(x_centers):
        row = free_energy[i]
        if len(row) != len(y_centers):
            raise ValueError("Each free_energy row must match y_centers length.")
        for j, y_val in enumerate(y_centers):
            rows.append((x_val, y_val, row[j]))
    _write_rows(output, ["x_center", "y_center", "free_energy"], rows)


def _write_pca(data, output):
    comps = data.get("principal_components", [])
    if not comps:
        _write_rows(output, ["frame"], [])
    else:
        width = len(comps[0])
        for row in comps:
            if len(row) != width:
                raise ValueError("principal_components rows must have the same width.")
        header = ["frame"] + [f"pc{i + 1}" for i in range(width)]
        rows = [(idx, *row) for idx, row in enumerate(comps)]
        _write_rows(output, header, rows)
    variance = data.get("variance")
    if isinstance(variance, list):
        variance_output = output[:-4] + "_variance.csv" if output.endswith(".csv") else output + "_variance.csv"
        _write_rows(variance_output, ["component", "cumulated_variance"],
                    [(idx + 1, value) for idx, value in enumerate(variance)])


def convert_json_to_csv(input_path, output_path):
    with open(input_path, "r") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("Only JSON objects are supported.")

    if "x" in data and "y" in data:
        _write_xy_series(data, output_path)
        return [output_path]
    if "times" in data and "files" in data:
        _write_times_files(data, output_path)
        return [output_path]
    if "free_energy" in data and "x_centers" in data and "y_centers" in data:
        _write_fes(data, output_path)
        return [output_path]
    if "principal_components" in data:
        _write_pca(data, output_path)
        outputs = [output_path]
        variance_output = output_path[:-4] + "_variance.csv" if output_path.endswith(".csv") else output_path + "_variance.csv"
        if isinstance(data.get("variance"), list):
            outputs.append(variance_output)
        return outputs

    raise ValueError(f"Unsupported JSON schema in {input_path}.")


def _iter_json_files(path, recursive):
    if os.path.isfile(path):
        return [path]
    results = []
    for root, _, files in os.walk(path):
        for name in files:
            if name.lower().endswith(".json"):
                results.append(os.path.join(root, name))
        if not recursive:
            break
    return sorted(results)


def json2csv(args):
    """
    This **function** converts supported Xponge JSON outputs to CSV files.

    :param args: arguments from argparse
    :return: None
    """
    input_path = args.input
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"can not find {input_path}")

    if os.path.isfile(input_path):
        output_path = args.output or os.path.splitext(input_path)[0] + ".csv"
        outputs = convert_json_to_csv(input_path, output_path)
        print("\n".join(outputs))
        return

    json_files = _iter_json_files(input_path, args.recursive)
    if not json_files:
        raise ValueError(f"can not find json files in {input_path}")

    outdir = args.output or input_path
    os.makedirs(outdir, exist_ok=True)
    generated = []
    for json_file in json_files:
        relpath = os.path.relpath(json_file, input_path)
        csv_path = os.path.join(outdir, os.path.splitext(relpath)[0] + ".csv")
        generated.extend(convert_json_to_csv(json_file, csv_path))
    print("\n".join(generated))
