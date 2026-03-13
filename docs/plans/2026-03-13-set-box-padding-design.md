# Design: Molecule.set_box_padding

Date: 2026-03-13

## Summary
Add a Molecule instance method `set_box_padding(padding=0.5, center=True)` to compute a tight periodic box from current coordinates, using a single user-defined padding as the only vacuum thickness. The method forces `box_length` to the computed values and optionally recenters coordinates so the molecule fits inside the box with the specified padding.

## Goals
- Provide an explicit API to make boxes compact by controlling padding only.
- Ensure coordinates are translated into the box when requested.
- Avoid implicit dependencies on `GlobalSetting.boxspace` or `vacuum_layer`.

## Non-Goals
- No changes to default box construction in `write_coordinate` or other IO paths.
- No automatic application during `load_pdb` or other loaders.
- No changes to `box_angle`.

## API
```python
Molecule.set_box_padding(padding=0.5, center=True)
```

### Behavior
1. Compute `min` and `max` over current atom coordinates.
2. Compute box lengths per axis: `(max - min) + 2 * padding`.
3. Force `self.box_length` to those values.
4. If `center=True`, translate all atom coordinates so the new minimum is exactly `padding` along each axis.

### Validation
- `padding < 0`: raise `ValueError`.
- No atoms present: raise `ValueError`.

### Compatibility
- `box_angle` remains unchanged.
- Does not consult `GlobalSetting.nocenter`, `GlobalSetting.boxspace`, or `vacuum_layer`.

## Tests
Add unit tests that:
- Verify `box_length == max-min + 2*padding` after calling the method.
- Verify min coordinate equals `padding` when `center=True`.
- Verify `padding < 0` raises `ValueError`.
- Verify empty molecule raises `ValueError`.

## Notes
This is an explicit opt-in API to let users tighten boxes for higher initial density without altering existing IO behavior.
