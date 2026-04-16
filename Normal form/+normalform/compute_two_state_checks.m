function checks = compute_two_state_checks(model, core, amTable, fmTable)
%COMPUTE_TWO_STATE_CHECKS Compute AM/FM residuals and structure-constraint checks.

checks = struct();

amResidual = normalform.frequency_from_linearization(model, amTable.f11, amTable.f12) - core.omega0;
fmResidual = normalform.first_lyapunov_from_tensors(model, fmTable.f11, fmTable.f12) - core.chi0;

checks.am_omega_max_abs = max(abs(amResidual), [], 'omitnan');
checks.fm_chi_max_abs = max(abs(fmResidual), [], 'omitnan');
checks.constraint_max_abs = max([ ...
    max(abs(amTable.f22 + amTable.f11), [], 'omitnan'), ...
    max(abs(fmTable.f22 + fmTable.f11), [], 'omitnan'), ...
    max(abs(amTable.f21 - core.a21_over_a12 .* amTable.f12), [], 'omitnan'), ...
    max(abs(fmTable.f21 - core.a21_over_a12 .* fmTable.f12), [], 'omitnan')], [], 'omitnan');
end
