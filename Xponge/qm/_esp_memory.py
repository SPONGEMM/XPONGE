"""Shared ESP memory-budget utilities."""

from __future__ import annotations

import math
import re


DEFAULT_ESP_MEMORY_LIMIT_BYTES = 1024 ** 3
DEFAULT_ESP_SAFETY_FACTOR = 0.8
DEFAULT_ESP_CHUNK_POLICY = "auto"
SUPPORTED_ESP_CHUNK_POLICIES = ("auto", "dual", "full", "grid", "pointwise")

_MEMORY_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([kmgt]?i?b?)?\s*$", re.IGNORECASE)
_MEMORY_SCALE = {
    "": 1,
    "b": 1,
    "k": 1024,
    "kb": 1024,
    "ki": 1024,
    "kib": 1024,
    "m": 1024 ** 2,
    "mb": 1024 ** 2,
    "mi": 1024 ** 2,
    "mib": 1024 ** 2,
    "g": 1024 ** 3,
    "gb": 1024 ** 3,
    "gi": 1024 ** 3,
    "gib": 1024 ** 3,
    "t": 1024 ** 4,
    "tb": 1024 ** 4,
    "ti": 1024 ** 4,
    "tib": 1024 ** 4,
}


def normalize_chunk_policy(policy):
    if policy is None:
        return DEFAULT_ESP_CHUNK_POLICY
    normalized = str(policy).strip().lower().replace("-", "_")
    aliases = {
        "grid_chunk": "grid",
        "shell_grid_chunk": "dual",
        "shellgridchunk": "dual",
        "ao_grid_chunk": "dual",
        "ao+grid": "dual",
        "point": "pointwise",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in SUPPORTED_ESP_CHUNK_POLICIES:
        supported = ", ".join(SUPPORTED_ESP_CHUNK_POLICIES)
        raise ValueError(f"ESP chunk policy should be one of: {supported}")
    return normalized


def parse_memory_limit_bytes(limit, default=DEFAULT_ESP_MEMORY_LIMIT_BYTES):
    if limit is None:
        return int(default)
    if isinstance(limit, (int, float)):
        parsed = int(limit)
        if parsed <= 0:
            raise ValueError("ESP memory limit should be positive")
        return parsed
    match = _MEMORY_RE.match(str(limit))
    if not match:
        raise ValueError(f"Invalid ESP memory limit: {limit!r}")
    value = float(match.group(1))
    unit = (match.group(2) or "").lower()
    scale = _MEMORY_SCALE.get(unit)
    if scale is None:
        raise ValueError(f"Invalid ESP memory limit unit: {limit!r}")
    parsed = int(value * scale)
    if parsed <= 0:
        raise ValueError("ESP memory limit should be positive")
    return parsed


def normalize_safety_factor(safety_factor, default=DEFAULT_ESP_SAFETY_FACTOR):
    if safety_factor is None:
        return float(default)
    factor = float(safety_factor)
    if not 0 < factor <= 1:
        raise ValueError("ESP safety factor should be in the range (0, 1]")
    return factor


def estimate_aux_tensor_bytes(naoi, naoj, naux, itemsize=8):
    return int(max(1, naoi) * max(1, naoj) * max(1, naux) * max(1, itemsize))


def iter_chunk_slices(length, chunk_size):
    chunk_size = max(1, int(chunk_size))
    for start in range(0, int(length), chunk_size):
        stop = min(int(length), start + chunk_size)
        yield start, stop


def build_shell_blocks(ao_loc, max_ao_per_block):
    ao_loc = list(ao_loc)
    max_ao_per_block = max(1, int(max_ao_per_block))
    blocks = []
    nbas = len(ao_loc) - 1
    shell_start = 0
    while shell_start < nbas:
        shell_stop = shell_start + 1
        while shell_stop < nbas and ao_loc[shell_stop + 1] - ao_loc[shell_start] <= max_ao_per_block:
            shell_stop += 1
        blocks.append((shell_start, shell_stop, ao_loc[shell_start], ao_loc[shell_stop]))
        shell_start = shell_stop
    return blocks


def choose_dual_chunk_layout(ao_loc, grid_count, usable_bytes, preferred_grid_chunk=128):
    grid_count = max(1, int(grid_count))
    usable_bytes = max(1, int(usable_bytes))
    preferred_grid_chunk = max(1, min(grid_count, int(preferred_grid_chunk)))
    max_ao_per_block = max(1, int(math.sqrt(usable_bytes / (8 * preferred_grid_chunk))))
    shell_blocks = build_shell_blocks(ao_loc, max_ao_per_block)
    largest_block_ao = max(block[3] - block[2] for block in shell_blocks)
    grid_chunk_size = max(1, min(grid_count, usable_bytes // max(1, 8 * largest_block_ao * largest_block_ao)))
    return shell_blocks, grid_chunk_size, largest_block_ao
