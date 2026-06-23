"""
This **module** is used to calculate the RESP charge.
The **module** is not available on Windows unless a supported backend is installed.
"""

from __future__ import annotations

from ..helper import Xprint, set_global_alternative_names
from ..qm import scheduler as qm_scheduler
from . import resp_core


def _normalize_backend_name(backend):
    try:
        return qm_scheduler.normalize_backend_name(backend)
    except ValueError as exc:
        raise ValueError(str(exc).replace("QM backend", "RESP backend", 1)) from exc


def _normalize_core_name(core):
    if core is None:
        return "python"
    core_name = str(core).strip().lower()
    if core_name != "python":
        raise ValueError("RESP core should be one of: python")
    return core_name


def resp_fit(assign, basis="6-31g*", opt=False, charge=None, spin=0, extra_equivalence=None,
             grid_density=6, grid_cell_layer=4, radius=None, a1=0.0005, a2=0.001,
             two_stage=True, only_esp=False, backend=None, core=None,
             esp_memory_limit=None, esp_chunk_policy="auto", esp_safety_factor=0.8):
    """
    This **function** fits the RESP partial charge for an Assign instance.

    :param assign: the Assign instance
    :param basis: the basis for Hartree-Fock calculation
    :param opt: whether do the geometry optimization
    :param charge: total charge of the molecule. If None, it will use the sum of the assign.charge
    :param spin: total spin of the molecule. ``S`` instead of ``2S+1``.
    :param extra_equivalence: the extra equivalence to constrain the charge
    :param grid_density: the density for grids to fit, in the unit of amgstrom^-3
    :param grid_cell_layer: the cell layer for grids to fit
    :param radius: the vdw radius for different elements. Default is ``default_radius`` in this module.
    :param a1: the restrain factor in the first step
    :param a2: the restrain factor in the second step
    :param two_stage: whether do the second stage fitting. If set to False, the second stage fitting will not be done
    :param only_esp: whether do the first stage fitting. If set to True, no restrained fitting will be done
    :param backend: the QM backend to use. Default is ``pyscf``.
    :param core: the RESP numerical core to use. Only ``python`` is currently supported in Xponge-origin.
    :return: a list of charges
    """
    if extra_equivalence is None:
        extra_equivalence = []
    if charge is None:
        charge = int(round(resp_core.np.sum(assign.charge)))

    backend_name = _normalize_backend_name(backend)
    _normalize_core_name(core)

    scf_result = qm_scheduler.run_scf(
        assign,
        backend=backend_name,
        basis=basis,
        charge=charge,
        spin=spin,
        optimize_geometry=opt,
    )
    grids = resp_core.get_mk_grid(
        assign,
        scf_result.coordinates_bohr,
        area_density=grid_density,
        layer=grid_cell_layer,
        radius=radius,
    )
    esp_result = qm_scheduler.compute_esp_on_grid(
        scf_result,
        grids,
        memory_limit=esp_memory_limit,
        chunk_policy=esp_chunk_policy,
        safety_factor=esp_safety_factor,
    )
    return resp_core.fit_resp_from_esp(
        assign,
        atom_coordinates_bohr=scf_result.coordinates_bohr,
        nuclear_charges=scf_result.nuclear_charges,
        grid_points_bohr=grids,
        esp_values_au=esp_result.electronic_esp_au,
        charge=charge,
        extra_equivalence=extra_equivalence,
        a1=a1,
        a2=a2,
        two_stage=two_stage,
        only_esp=only_esp,
    )


_get_mk_grid = resp_core._get_mk_grid
_resp_scf_kernel = resp_core._resp_scf_kernel
_find_tofit_second = resp_core._find_tofit_second
_correct_extra_equivalence = resp_core._correct_extra_equivalence
_get_a20_and_b20 = resp_core._get_a20_and_b20
_find_restrained_second_stage_groups = resp_core._find_restrained_second_stage_groups

set_global_alternative_names()

Xprint("""Reference for resp.py:
1. pyscf
  Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu, J. McClain, S. Sharma, S. Wouters, and G. K.-L. Chan
    PySCF: the Python-based simulations of chemistry framework
    WIREs Computational Molecular Science 2018 8(e1340)
    DOI: 10.1002/wcms.1340

2. ESP MK grid generation
  Brent H. Besler, Kenneth M. Merz Jr., Peter A. Kollman
    Atomic charges derived from semiempirical methods
    Journal of Computational Chemistry 1990 11 431-439
    DOI: 10.1002/jcc.540110404

3. RESP
  Christopher I. Bayly, Piotr Cieplak, Wendy Cornell, and Peter A. Kollman
    A well-behaved electrostatic potential-based method using charge restraints for deriving atomic char
    Journal of Physical Chemistry 1993 97(40) 10269-10280
    DOI: 10.1021/j100142a004

""")
