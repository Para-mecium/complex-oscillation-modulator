function controlTable = make_control_table(f11, f12, a21_over_a12)
%MAKE_CONTROL_TABLE Build a standard control-curve table for two-state systems.

controls = normalform.complete_two_state_controls(f11, f12, a21_over_a12);
controlTable = array2table(controls, ...
    'VariableNames', {'f11', 'f12', 'f21', 'f22'});
end
