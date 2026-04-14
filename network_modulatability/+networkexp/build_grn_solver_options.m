function solverOptions = build_grn_solver_options(params)
%BUILD_GRN_SOLVER_OPTIONS Assemble sparse solver metadata for GRN networks.

d_p = params{1};
d_E = params{2};
k_m = params{3};
epsilon = params{4};
S = params{5};
K_d = params{6};
n = params{7};
d_m = params{8};
matAdj = sparse(params{9});

numNodes = size(matAdj, 1);
identity = speye(numNodes);
jPattern = [spones(matAdj) + identity, identity; identity, identity];

solverOptions = struct();
solverOptions.JPattern = jPattern;
solverOptions.Jacobian = @(t, x) evaluate_grn_jacobian(t, x, d_p, d_E, k_m, epsilon, S, K_d, n, d_m, matAdj); %#ok<INUSD>
end

function jacobian = evaluate_grn_jacobian(~, y, d_p, d_E, k_m, epsilon, S, K_d, n, d_m, matAdj)
y = y(:);
numNodes = size(matAdj, 1);
p = y(1:numNodes);
m = y(numNodes + 1:end);

den = k_m + p + p.^2;
g = p ./ (1 + p);
gPrime = 1 ./ (1 + p).^2;
hillDen = K_d^n + p.^n;
hillPrime = -S * K_d^n * n * p.^(n - 1) ./ (hillDen.^2);

jpp = spdiags(-d_p - d_E * (k_m - p.^2) ./ (den.^2), 0, numNodes, numNodes) + ...
    spdiags(m, 0, numNodes, numNodes) * matAdj * spdiags(gPrime, 0, numNodes, numNodes);
jpm = spdiags(1 + matAdj * g, 0, numNodes, numNodes);
jmp = spdiags(epsilon * hillPrime, 0, numNodes, numNodes);
jmm = -epsilon * d_m * speye(numNodes);

jacobian = [jpp, jpm; jmp, jmm];
end
