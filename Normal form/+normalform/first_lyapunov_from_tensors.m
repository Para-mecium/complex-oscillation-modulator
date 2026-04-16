function chi = first_lyapunov_from_tensors(model, f11, f12)
%FIRST_LYAPUNOV_FROM_TENSORS Compute chi from A/H/T tensors for a two-state system.

chi = evaluate_two_input_scalar(@(x, y) chi_scalar(model, x, y), f11, f12);
end

function chi = chi_scalar(model, f11, f12)
A = model.A_function(f11, f12);
H = model.H_function(f11, f12);
T = model.T_function(f11, f12);

try
    [omega, q0, p0] = matcont_hopf_basis(A);
    B = @(x, y) apply_quadratic_tensor(H, x, y);
    C = @(x, y, z) apply_cubic_tensor(T, x, y, z);

    h20 = (2i * omega * eye(2) - A) \ B(q0, q0);
    h11 = -A \ B(q0, conj(q0));
    h21 = C(q0, q0, conj(q0)) + 2 * B(q0, h11) + B(h20, conj(q0));
    chi = real(p0' * h21) / 2;
catch err
    if strcmp(err.identifier, 'normalform:NoHopfPair')
        chi = nan;
        return
    end
    rethrow(err);
end
end

function [omega, q0, p0] = matcont_hopf_basis(A)
[rightVectors, rightValues] = eig(A);
[idx1, idx2] = select_hopf_pair(diag(rightValues));
[basisV, ~] = qr([real(rightVectors(:, idx1)), imag(rightVectors(:, idx1))], 0);

[leftVectors, leftValues] = eig(A');
[leftIdx1, ~] = select_hopf_pair(diag(leftValues));
[basisW, ~] = qr([real(leftVectors(:, leftIdx1)), imag(leftVectors(:, leftIdx1))], 0);

k = real(rightValues(idx1, idx1) * rightValues(idx2, idx2));
if k <= 0
    error('normalform:NoHopfPair', ...
        'Controlled linearization does not define a positive Hopf frequency.');
end

omega = sqrt(k);
reduced = A * A + k * eye(size(A));
bordered = [reduced, basisW; basisV', zeros(2)];
rhs = [zeros(size(A, 1), 2); eye(2)];
vext = bordered \ rhs;
wext = bordered' \ rhs;

alpha = vext(1:2, 1)' * A * vext(1:2, 2) - 1i * omega * (vext(1:2, 1)' * vext(1:2, 2));
beta = -vext(1:2, 1)' * A * vext(1:2, 1) + 1i * omega * (vext(1:2, 1)' * vext(1:2, 1));
q0 = alpha * vext(1:2, 1) + beta * vext(1:2, 2);

alpha = wext(1:2, 1)' * A' * wext(1:2, 2) + 1i * omega * (wext(1:2, 1)' * wext(1:2, 2));
beta = -wext(1:2, 1)' * A' * wext(1:2, 1) - 1i * omega * (wext(1:2, 1)' * wext(1:2, 1));
p0 = alpha * wext(1:2, 1) + beta * wext(1:2, 2);

q0 = q0 / norm(q0);
p0 = p0 / (q0' * p0);
end

function [idx1, idx2] = select_hopf_pair(eigenvalues)
smallestSum = inf;
idx1 = 1;
idx2 = 2;
for i = 1:(numel(eigenvalues) - 1)
    [candidate, offset] = min(abs(eigenvalues((i + 1):end) + eigenvalues(i)));
    if candidate < smallestSum
        smallestSum = candidate;
        idx1 = i;
        idx2 = i + offset;
    end
end
end

function value = apply_quadratic_tensor(H, x, y)
value = zeros(2, 1);
for i = 1:2
    for j = 1:2
        for k = 1:2
            value(i) = value(i) + H(i, j, k) * x(j) * y(k);
        end
    end
end
end

function value = apply_cubic_tensor(T, x, y, z)
value = zeros(2, 1);
for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                value(i) = value(i) + T(i, j, k, l) * x(j) * y(k) * z(l);
            end
        end
    end
end
end

function values = evaluate_two_input_scalar(fn, f11, f12)
if ~isequal(size(f11), size(f12))
    error('normalform:InputSizeMismatch', ...
        'f11 and f12 must have the same size.');
end

values = arrayfun(fn, f11, f12);
end
