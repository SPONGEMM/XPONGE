"""
Post-analysis tools for SPONGE trajectories.
"""
import hashlib
import inspect
import json
import os
import shlex

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_Sponge_trajectory(topo, traj, box):
    """
    Load a Sponge trajectory using MDAnalysis.

    :param topo: Path to the topology file.
    :param traj: Path to the trajectory file.
    :param box: Path to the box file.
    :return: MDAnalysis Universe object.
    """
    import MDAnalysis as mda
    import Xponge.analysis.md_analysis as xmda

    traj_list = _split_paths(traj)
    if not traj_list:
        raise ValueError("Trajectory file is required.")

    if len(traj_list) == 1:
        traj_path = traj_list[0]
        if box:
            if isinstance(box, (list, tuple)) and len(box) == 1:
                box = box[0]
            return mda.Universe(topo, traj_path, format=xmda.SpongeTrajectoryReader, box=box)
        return mda.Universe(topo, traj_path, format=xmda.SpongeTrajectoryReader)

    if not box:
        raise ValueError("Box file is required for Sponge trajectories.")
    if isinstance(box, (list, tuple)):
        if len(box) == 1:
            return mda.Universe(topo, traj_list, format=xmda.SpongeTrajectoryReader, box=box[0])
        if len(box) != len(traj_list):
            raise ValueError("Number of box files must match trajectories or be a single box file.")
        filenames = [
            (traj, xmda.SpongeTrajectoryReader.with_arguments(box=box_item))
            for traj, box_item in zip(traj_list, box)
        ]
        return mda.Universe(topo, filenames)
    return mda.Universe(topo, traj_list, format=xmda.SpongeTrajectoryReader, box=box)


def _selection_suffix(selection: str) -> str:
    value = (selection or "all").strip() or "all"
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:6]


def _traj_label(traj_path: str, index: int) -> str:
    base = os.path.basename(traj_path) if traj_path else ""
    name, _ = os.path.splitext(base)
    name = name.strip() or f"traj{index + 1}"
    if traj_path:
        suffix = hashlib.sha1(traj_path.encode("utf-8")).hexdigest()[:6]
    else:
        suffix = hashlib.sha1(str(index).encode("utf-8")).hexdigest()[:6]
    return f"{name}-{suffix}"


def _split_paths(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    if isinstance(value, str):
        if "," in value:
            return [item.strip() for item in value.split(",") if item.strip()]
        return [value.strip()] if value.strip() else []
    return [str(value).strip()] if str(value).strip() else []


def json_filename_for_selection(base_name: str, selection: str) -> str:
    return f"{base_name}_{_selection_suffix(selection)}.json"


def plot_filename_for_selection(base_name: str, selection: str) -> str:
    return f"{base_name}_{_selection_suffix(selection)}.png"


def _rmsd_values(rmsd_analysis):
    if hasattr(rmsd_analysis, "results") and hasattr(rmsd_analysis.results, "rmsd"):
        return rmsd_analysis.results.rmsd
    return rmsd_analysis.rmsd


def perform_pca(u, n_components=2, selection="backbone", outdir="."):
    """
    Perform PCA, save component plots, and print JSON.
    """
    from MDAnalysis.analysis.pca import PCA

    pca = PCA(u, select=selection, n_components=n_components)
    pca.run()
    comps = pca.results.p_components
    var = pca.results.cumulated_variance
    for i in range(comps.shape[1]):
        plt.figure()
        plt.plot(comps[:, i])
        plt.xlabel("Frame")
        plt.ylabel(f"PC{i+1}")
        plt.title(f"PCA Component {i+1}")
        plt.savefig(os.path.join(outdir, f"pca_pc{i+1}.png"))
        plt.close()
    result = {"principal_components": comps.tolist(), "variance": var.tolist()}
    with open(os.path.join(outdir, "pca_result.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result


def compute_free_energy_surface(u, cv1_func, cv2_func, bins=50, temperature=300, outdir="."):
    kB = 0.0019872041
    data1 = cv1_func(u)
    data2 = cv2_func(u)
    H, x_edges, y_edges = np.histogram2d(data1, data2, bins=bins, density=True)
    fe = -kB * temperature * np.log(H + 1e-12)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    plt.figure()
    plt.contourf(x_centers, y_centers, fe.T)
    plt.xlabel("CV1")
    plt.ylabel("CV2")
    plt.title("Free Energy Surface")
    plt.colorbar(label="Free Energy (kcal/mol)")
    plt.savefig(os.path.join(outdir, "free_energy_surface.png"))
    plt.close()
    result = {
        "free_energy": fe.tolist(),
        "x_centers": x_centers.tolist(),
        "y_centers": y_centers.tolist(),
    }
    with open(os.path.join(outdir, "free_energy_surface.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result


def compute_rmsf(u, outdir=".", selection="backbone and name CA"):
    from MDAnalysis.analysis.rms import RMSF

    os.makedirs(outdir, exist_ok=True)
    ag = u.select_atoms(selection)
    rmsf_analysis = RMSF(ag)
    rmsf_analysis.run()
    values = rmsf_analysis.results.rmsf
    plt.figure()
    plt.plot(values)
    plt.xlabel("Atom Index")
    plt.ylabel("RMSF (Å)")
    plt.title(f"RMSF of {selection}")
    plt.savefig(os.path.join(outdir, plot_filename_for_selection("rmsf", selection)))
    plt.close()
    x = np.arange(len(values))
    result = {"selection": selection, "x": x.tolist(), "y": values.tolist()}
    with open(os.path.join(outdir, json_filename_for_selection("rmsf", selection)), "w") as f:
        json.dump(result, f, indent=2)
    return result


def rmsd_calculator(topo_path, traj_path, box_path, selection, delta_t_ns, outdir=".", universe=None):
    from MDAnalysis.analysis import rms

    os.makedirs(outdir, exist_ok=True)
    u = universe or load_Sponge_trajectory(topo_path, traj_path, box_path)
    ag = u.select_atoms(selection)
    rmsd = rms.RMSD(ag, superposition=True, ref_frame=0)
    rmsd.run()
    rmsds = np.array(_rmsd_values(rmsd)[:, 2])
    plt.figure()
    time = np.arange(len(rmsds)) * delta_t_ns
    plt.plot(time, rmsds, label="RMSD")
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD over Time")
    plt.savefig(os.path.join(outdir, plot_filename_for_selection("rmsd", selection)))
    plt.close()
    total_step = len(rmsds)
    last_quarter = int(3 * total_step / 4)
    tail = rmsds[last_quarter:] if total_step else np.array([])
    stats = {
        "avg": float(rmsds.mean()) if total_step else 0.0,
        "last_quarter_mean": float(tail.mean()) if tail.size else 0.0,
        "last_quarter_std": float(tail.std()) if tail.size else 0.0,
        "total_frames": int(total_step),
    }
    result = {"selection": selection, "x": time.tolist(), "y": rmsds.tolist(), "stats": stats}
    with open(os.path.join(outdir, json_filename_for_selection("rmsd", selection)), "w") as f:
        json.dump(result, f, indent=2)
    return result


def extract_pdb(topo_path, traj_path, box_path, times, delta_t_ns, selection="all", outdir=".", universe=None):
    os.makedirs(outdir, exist_ok=True)
    u = universe or load_Sponge_trajectory(topo_path, traj_path, box_path)
    n_frames = len(u.trajectory)
    ag = u.select_atoms(selection)
    files = []
    times_used = []
    for t in times:
        step = int(round(t / delta_t_ns))
        if step < 0 or step >= n_frames:
            continue
        u.trajectory[step]
        file = os.path.join(outdir, f"time_{t}ns.pdb")
        files.append(file)
        times_used.append(t)
        ag.write(file)
    result = {"selection": selection, "files": files, "times": times_used}
    with open(os.path.join(outdir, json_filename_for_selection("extract_pdb", selection)), "w") as f:
        json.dump(result, f, indent=2)
    return result


def compute_radius_of_gyration(u, delta_t_ps, outdir=".", selection="protein and backbone"):
    os.makedirs(outdir, exist_ok=True)
    ag = u.select_atoms(selection)
    u.trajectory[0]
    rg = [ag.radius_of_gyration() for _ in u.trajectory]
    plt.figure()
    delta_t_ns = delta_t_ps / 1000.0
    time = np.arange(len(rg)) * delta_t_ns
    plt.plot(time, rg, label="Radius of Gyration")
    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Radius of Gyration vs Time")
    plt.savefig(os.path.join(outdir, plot_filename_for_selection("radius_of_gyration", selection)))
    plt.close()
    result = {"selection": selection, "x": time.tolist(), "y": [float(val) for val in rg]}
    with open(os.path.join(outdir, json_filename_for_selection("radius_of_gyration", selection)), "w") as f:
        json.dump(result, f, indent=2)
    return result


def analyze_hydrogen_bonds(u, delta_t_ps, between, selection="all", outdir=".", update_select=None):
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

    os.makedirs(outdir, exist_ok=True)
    default_between = ["protein", "resname SOL"]
    if between is None:
        between = []
    if isinstance(between, str):
        between = [entry.strip() for entry in between.split(",") if entry.strip()]
    between = [str(entry).strip() for entry in (between or []) if str(entry).strip()]
    if not between:
        between = default_between

    if update_select is None:
        update_select_flag = False if between == default_between else True
    else:
        update_select_flag = bool(update_select)
    kwargs = {"universe": u, "between": between}
    try:
        signature = inspect.signature(HydrogenBondAnalysis)
        if "update_selections" in signature.parameters:
            kwargs["update_selections"] = update_select_flag
        elif "update_selection" in signature.parameters:
            kwargs["update_selection"] = update_select_flag
    except Exception:
        kwargs["update_selections"] = update_select_flag

    hba = HydrogenBondAnalysis(**kwargs)
    hba.run()
    hbonds = hba.results.hbonds.tolist()
    frames = np.arange(u.trajectory.n_frames)
    delta_t_ns = delta_t_ps / 1000.0
    time = frames * delta_t_ns
    counts = hba.count_by_time()
    plt.plot(time, counts)
    plt.xlabel("Time (ns)")
    plt.ylabel("Number of Hydrogen Bonds")
    plt.title("Hydrogen Bonds Over Time")
    plt.savefig(os.path.join(outdir, plot_filename_for_selection("hydrogen_bonds", selection)))
    plt.close()
    result = {"selection": selection, "x": time.tolist(), "y": counts.tolist(), "hbonds": hbonds}
    with open(os.path.join(outdir, json_filename_for_selection("hydrogen_bonds", selection)), "w") as f:
        json.dump(result, f, indent=2)
    return result


def _resolve_box_path(traj_path, box_path):
    traj_list = _split_paths(traj_path)
    if not traj_list:
        return box_path

    box_list = _split_paths(box_path) if box_path is not None else []
    if box_list:
        if len(box_list) == 1 or len(box_list) == len(traj_list):
            return box_list[0] if len(box_list) == 1 else box_list
        raise ValueError("Number of box files must match trajectories or be a single box file.")

    resolved = []
    for traj in traj_list:
        dirname, basename = os.path.split(traj)
        if basename == "mdcrd.dat":
            candidate = "mdbox.txt"
        else:
            candidate = basename.replace(".dat", ".box")
        candidate = os.path.join(dirname, candidate)
        resolved.append(candidate if os.path.exists(candidate) else None)

    if len(traj_list) == 1:
        return resolved[0]
    missing = [traj_list[i] for i, item in enumerate(resolved) if item is None]
    if missing:
        missing_list = ", ".join(missing)
        raise ValueError(f"Box file not found for trajectories: {missing_list}. Please provide -b/--box.")
    return resolved


def _resolve_dt_ns(dt_ns, dt_ps):
    if dt_ns is not None and dt_ps is not None:
        raise ValueError("Please set only one of --dt-ns or --dt-ps.")
    if dt_ns is not None:
        return dt_ns
    if dt_ps is not None:
        return dt_ps / 1000.0
    raise ValueError("Please provide --dt-ns or --dt-ps.")


def _resolve_dt_ps(dt_ps, dt_ns):
    if dt_ns is not None and dt_ps is not None:
        raise ValueError("Please set only one of --dt-ps or --dt-ns.")
    if dt_ps is not None:
        return dt_ps
    if dt_ns is not None:
        return dt_ns * 1000.0
    raise ValueError("Please provide --dt-ps or --dt-ns.")


def _cv_from_spec(spec):
    if not spec:
        raise ValueError("CV spec is empty.")
    name, _, arg = spec.partition(":")
    name = name.strip().lower()
    selection = arg.strip() if arg else ""

    if name == "rmsd":
        if not selection:
            selection = "backbone"

        def _cv(u):
            from MDAnalysis.analysis import rms

            ag = u.select_atoms(selection)
            rmsd = rms.RMSD(ag, superposition=True, ref_frame=0)
            rmsd.run()
            return np.array(_rmsd_values(rmsd)[:, 2])

        return _cv
    if name == "rgyr":
        if not selection:
            selection = "protein and backbone"

        def _cv(u):
            ag = u.select_atoms(selection)
            u.trajectory[0]
            return np.array([ag.radius_of_gyration() for _ in u.trajectory])

        return _cv
    raise ValueError(f"Unsupported CV spec: {spec}. Use rmsd:SELECTION or rgyr:SELECTION.")


def _parse_times(value):
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    if value is None:
        return []
    return [float(v) for v in str(value).split(",") if str(v).strip()]


def _run_script(u, args):
    with open(args.input, "r") as f:
        lines = f.readlines()
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = shlex.split(line)
        cmd = tokens[0].lower()
        kwargs = {}
        for token in tokens[1:]:
            if "=" not in token:
                raise ValueError(f"Invalid token '{token}' in script line: {raw_line}")
            key, value = token.split("=", 1)
            kwargs[key.strip().lower()] = value
        outdir = kwargs.get("outdir", args.outdir)
        traj_spec = args.traj
        box_spec = args.box
        line_traj = kwargs.get("traj")
        line_box = kwargs.get("box")
        if line_traj is not None or line_box is not None:
            traj_spec = line_traj if line_traj is not None else args.traj
            box_spec = line_box if line_box is not None else args.box
            box_spec = _resolve_box_path(traj_spec, box_spec)
            line_u = load_Sponge_trajectory(args.topo, traj_spec, box_spec)
        else:
            line_u = u
        if cmd == "rmsd":
            selection = kwargs.get("selection", "backbone")
            dt_ns = float(kwargs["dt_ns"]) if "dt_ns" in kwargs else args.dt_ns
            dt_ps = float(kwargs["dt_ps"]) if "dt_ps" in kwargs else args.dt_ps
            delta_t_ns = _resolve_dt_ns(dt_ns, dt_ps)
            line_u.trajectory[0]
            rmsd_calculator(args.topo, traj_spec, box_spec, selection, delta_t_ns, outdir=outdir, universe=line_u)
        elif cmd == "rmsf":
            selection = kwargs.get("selection", "backbone and name CA")
            line_u.trajectory[0]
            compute_rmsf(line_u, outdir=outdir, selection=selection)
        elif cmd == "rgyr":
            selection = kwargs.get("selection", "protein and backbone")
            dt_ps = float(kwargs["dt_ps"]) if "dt_ps" in kwargs else args.dt_ps
            dt_ns = float(kwargs["dt_ns"]) if "dt_ns" in kwargs else args.dt_ns
            delta_t_ps = _resolve_dt_ps(dt_ps, dt_ns)
            line_u.trajectory[0]
            compute_radius_of_gyration(line_u, delta_t_ps, outdir=outdir, selection=selection)
        elif cmd == "hbond":
            selection = kwargs.get("selection", "all")
            between = kwargs.get("between", None)
            dt_ps = float(kwargs["dt_ps"]) if "dt_ps" in kwargs else args.dt_ps
            dt_ns = float(kwargs["dt_ns"]) if "dt_ns" in kwargs else args.dt_ns
            delta_t_ps = _resolve_dt_ps(dt_ps, dt_ns)
            update_select = kwargs.get("update_select", None)
            line_u.trajectory[0]
            analyze_hydrogen_bonds(
                line_u,
                delta_t_ps,
                between,
                selection=selection,
                outdir=outdir,
                update_select=update_select,
            )
        elif cmd == "pca":
            n_components = int(kwargs.get("n_components", kwargs.get("n", 2)))
            selection = kwargs.get("selection", "backbone")
            line_u.trajectory[0]
            perform_pca(line_u, n_components=n_components, selection=selection, outdir=outdir)
        elif cmd == "fes":
            cv1 = kwargs.get("cv1", "")
            cv2 = kwargs.get("cv2", "")
            bins = int(kwargs.get("bins", 50))
            temp = float(kwargs.get("temperature", kwargs.get("temp", 300)))
            line_u.trajectory[0]
            compute_free_energy_surface(
                line_u,
                _cv_from_spec(cv1),
                _cv_from_spec(cv2),
                bins=bins,
                temperature=temp,
                outdir=outdir,
            )
        elif cmd == "extract_pdb":
            selection = kwargs.get("selection", "all")
            times = _parse_times(kwargs.get("times", ""))
            dt_ns = float(kwargs["dt_ns"]) if "dt_ns" in kwargs else args.dt_ns
            dt_ps = float(kwargs["dt_ps"]) if "dt_ps" in kwargs else args.dt_ps
            delta_t_ns = _resolve_dt_ns(dt_ns, dt_ps)
            line_u.trajectory[0]
            extract_pdb(
                args.topo,
                traj_spec,
                box_spec,
                times,
                delta_t_ns,
                selection=selection,
                outdir=outdir,
                universe=line_u,
            )
        else:
            raise ValueError(f"Unsupported command in script: {cmd}")


def run_traj_cli(args):
    args.traj = _split_paths(args.traj)
    box_path = _resolve_box_path(args.traj, args.box)
    args.box = box_path
    os.makedirs(args.outdir, exist_ok=True)
    traj_list = _split_paths(args.traj)
    if not traj_list:
        raise ValueError("Trajectory file is required.")
    if isinstance(args.box, (list, tuple)):
        box_list = list(args.box)
    else:
        box_list = [args.box] * len(traj_list)

    def _run_single(traj_path, box_path, outdir):
        u = load_Sponge_trajectory(args.topo, traj_path, box_path)
        if args.input:
            _run_script(u, args)
            return
        if not args.analysis_cmd:
            raise ValueError("Please provide a subcommand or use -i/--input.")
        if args.analysis_cmd == "rmsd":
            delta_t_ns = _resolve_dt_ns(args.dt_ns, args.dt_ps)
            u.trajectory[0]
            rmsd_calculator(args.topo, traj_path, box_path, args.selection, delta_t_ns, outdir=outdir, universe=u)
        elif args.analysis_cmd == "rmsf":
            u.trajectory[0]
            compute_rmsf(u, outdir=outdir, selection=args.selection)
        elif args.analysis_cmd == "rgyr":
            delta_t_ps = _resolve_dt_ps(args.dt_ps, args.dt_ns)
            u.trajectory[0]
            compute_radius_of_gyration(u, delta_t_ps, outdir=outdir, selection=args.selection)
        elif args.analysis_cmd == "hbond":
            delta_t_ps = _resolve_dt_ps(args.dt_ps, args.dt_ns)
            u.trajectory[0]
            analyze_hydrogen_bonds(
                u,
                delta_t_ps,
                args.between,
                selection=args.selection,
                outdir=outdir,
                update_select=args.update_select,
            )
        elif args.analysis_cmd == "pca":
            u.trajectory[0]
            perform_pca(u, n_components=args.n_components, selection=args.selection, outdir=outdir)
        elif args.analysis_cmd == "fes":
            u.trajectory[0]
            compute_free_energy_surface(
                u,
                _cv_from_spec(args.cv1),
                _cv_from_spec(args.cv2),
                bins=args.bins,
                temperature=args.temperature,
                outdir=outdir,
            )
        elif args.analysis_cmd == "extract_pdb":
            delta_t_ns = _resolve_dt_ns(args.dt_ns, args.dt_ps)
            u.trajectory[0]
            extract_pdb(
                args.topo,
                traj_path,
                box_path,
                args.times,
                delta_t_ns,
                selection=args.selection,
                outdir=outdir,
                universe=u,
            )
        else:
            raise ValueError(f"Unsupported subcommand: {args.analysis_cmd}")

    if len(traj_list) > 1 and not args.input:
        for idx, traj_path in enumerate(traj_list):
            box_path = box_list[idx] if idx < len(box_list) else None
            label = _traj_label(traj_path, idx)
            outdir = os.path.join(args.outdir, label)
            os.makedirs(outdir, exist_ok=True)
            _run_single(traj_path, box_path, outdir)
        return

    _run_single(traj_list[0], box_list[0], args.outdir)
