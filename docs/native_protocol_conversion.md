# Native protocol conversion

`convert_bundle_to_legacy` exports named positional restraints (including the
full-system reference coordinates), CV harmonic restraints and their schedules,
steering, and SITS method defaults. Disabled objects are omitted. Explicit mdin
settings retain priority over native defaults. SITS defaults are emitted as mdin
commands; selection and Nk are separate legacy files. Legacy SITS has no input
for checkpoint `log_norm` or `log_nk`, so exporting those states raises an error
rather than silently discarding them.

Legacy-to-bundle conversion no longer copies or binds text sidecars after
successful typed parsing of native-supported topology, CV/restraint/steering,
SITS, constraints, soft walls, and Nose–Hoover inputs. Dynamic custom-force data
also lives only in the H5 topology. Legacy CV-family configuration remains in
H5 `/config` datasets and is read directly by SPONGE. Reading `/steer/config`
and legacy `steer_cv_in_file` requires the accompanying SPONGE update.

Fallback sidecars remain for unsupported raw imports, QC type input and legacy
MetaD edge/checkpoint formats whose native runtime state is not yet equivalent.
Reverse export remains compatible with older bundles that contain sidecars.
