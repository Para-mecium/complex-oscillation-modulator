# `orbitmatcont`

`orbitmatcont.extract_periodic_orbit` is a new MATCONT-backed periodic-orbit interface.

It is intentionally separate from `PO_extract/extract_periodic_orbit.m`.

MATCONT-specific incremental defaults live in
[`+orbitmatcont/default_config.m`](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/+orbitmatcont/default_config.m).
Concrete task scripts should override only the fields they need.

## Call shape

```matlab
result = orbitmatcont.extract_periodic_orbit(odefun, y0, parameter, opts)
```

The four-argument signature is preserved, but the option contract is different from the
direct backend.

## What it does

1. forward integrate a pilot trajectory
2. estimate a rough period from the pilot tail
3. integrate a shorter orbit window of length `~ 1.25 * T_hat`
4. call MATCONT `initOrbLC` on that window
5. run a one-step `cont(@limitcycle, ...)` correction onto the LC branch
6. continue the LC branch with a MATCONT userfunction `p_active - p_input`
7. read the userfunction zero from MATCONT special points instead of guessing a branch column
8. accept the orbit only if the located active parameter matches the input value within tolerance

If automatic period estimation is not reliable for a given model, set:

- `opts.matcont_window_timespan`

This bypasses the automatic short-window choice and follows the documented
`initOrbLC` workflow directly.

## Minimum MATCONT-specific options

- `opts.matcont_active_parameter`
- `opts.matcont_odefile`

`opts.matcont_odefile` may be omitted only when the first argument `odefun` is itself a
MATCONT odefile.

If `parameter` is not the numeric parameter vector expected by MATCONT, provide:

- `opts.matcont_parameter_values`

## Output intent

The returned `result` keeps the same top-level stable fields used by the current project:

- `success`
- `has_orbit`
- `status`
- `code`
- `message`
- `period`
- `orbit_t`
- `orbit_y`
- `amplitude`
- `max_variable`
- `min_variable`

It also reports parameter-consistency fields:

- `input_parameter_values`
- `output_parameter_values`
- `active_parameter_index`
- `input_active_parameter_value`
- `seed_corrected_parameter_value`
- `output_active_parameter_value`
- `parameter_error`
- `parameter_status`

`orbitmatcont` treats parameter consistency as a hard requirement:

- success means the returned LC matches the input active parameter within `opts.matcont_parameter_tolerance`
- if MATCONT continues the branch but does not locate the userfunction zero, the interface throws `orbitmatcont:ParameterReturnFailed`
- MATCONT and solver failures are no longer wrapped into status codes; their native errors are allowed to surface

## Feasibility note

The current four-argument function signature can be preserved.

What cannot be preserved transparently is the old assumption that a generic ODE RHS alone
is enough. MATCONT still needs either:

- a MATCONT odefile, or
- a wrapper path that can generate one.

This implementation chooses the first path and keeps the extra MATCONT requirements inside
`opts`.
