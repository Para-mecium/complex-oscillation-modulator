function omega = frequency_from_linearization(model, f11, f12)
%FREQUENCY_FROM_LINEARIZATION Compute Hopf frequency from the controlled linearization.

omega = evaluate_two_input_scalar(@(x, y) frequency_scalar(model, x, y), f11, f12);
end

function omega = frequency_scalar(model, f11, f12)
A = model.A_function(f11, f12);
traceA = trace(A);
detA = det(A);
omegaSq = detA - 0.25 * traceA.^2;

if ~isfinite(omegaSq) || real(omegaSq) < 0
    omega = nan;
    return
end

omega = sqrt(real(omegaSq));
end

function values = evaluate_two_input_scalar(fn, f11, f12)
if ~isequal(size(f11), size(f12))
    error('normalform:InputSizeMismatch', ...
        'f11 and f12 must have the same size.');
end

values = arrayfun(fn, f11, f12);
end
