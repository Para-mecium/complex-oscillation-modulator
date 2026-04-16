function controls = complete_two_state_controls(f11, f12, a21_over_a12)
%COMPLETE_TWO_STATE_CONTROLS Recover f21 and f22 from the two structure constraints.

f11 = f11(:);
f12 = f12(:);

if numel(f11) ~= numel(f12)
    error('normalform:SizeMismatch', ...
        'f11 and f12 must contain the same number of points.');
end

controls = [f11, f12, a21_over_a12 .* f12, -f11];
end
