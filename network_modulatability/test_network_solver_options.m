function test_network_solver_options()
%TEST_NETWORK_SOLVER_OPTIONS Verify sparse Jacobian metadata for network models.

networkDir = fileparts(mfilename('fullpath'));
addpath(networkDir);

rng(1);
test_fhn_solver_options();
test_grn_solver_options();
test_nonstiff_solver_ignores_stiff_metadata();

disp('All network solver option tests passed.');
end

function test_fhn_solver_options()
numNodes = 4;
weightMatrix = sparse([0 1 0 0; 0 0 1 0; 1 0 0 1; 0 0 0 0]);
params = {0.5, 1, 1e-2, 1.3, 0.27, weightMatrix};
state = linspace(-0.2, 0.3, 2 * numNodes).';

solverOptions = networkexp.build_fhn_solver_options(params);
assert(isfield(solverOptions, 'JPattern') && isfield(solverOptions, 'Jacobian'), ...
    'FHN solver options must provide both JPattern and Jacobian.');

analyticJacobian = solverOptions.Jacobian(0, state);
finiteDifferenceJacobian = finite_difference(@(x) fhn_rhs(x, params), state);
assert(norm(analyticJacobian - finiteDifferenceJacobian, inf) < 1e-5, ...
    'FHN analytic Jacobian does not match finite differences.');
assert(all(all(spones(analyticJacobian) <= spones(solverOptions.JPattern))), ...
    'FHN JPattern must include every nonzero Jacobian entry.');
end

function test_grn_solver_options()
numNodes = 4;
weightMatrix = sparse([0 1 0 0; 0 0 1 0; 1 0 0 0; 0 1 0 0]);
params = {0.01, 1, 0.5, 0.05, 0.7, 1, 4, 1, weightMatrix};
state = linspace(0.2, 1.1, 2 * numNodes).';

solverOptions = networkexp.build_grn_solver_options(params);
assert(isfield(solverOptions, 'JPattern') && isfield(solverOptions, 'Jacobian'), ...
    'GRN solver options must provide both JPattern and Jacobian.');

analyticJacobian = solverOptions.Jacobian(0, state);
finiteDifferenceJacobian = finite_difference(@(x) grn_rhs(x, params), state);
assert(norm(analyticJacobian - finiteDifferenceJacobian, inf) < 1e-5, ...
    'GRN analytic Jacobian does not match finite differences.');
assert(all(all(spones(analyticJacobian) <= spones(solverOptions.JPattern))), ...
    'GRN JPattern must include every nonzero Jacobian entry.');
end

function test_nonstiff_solver_ignores_stiff_metadata()
solverOptions = struct( ...
    'JPattern', speye(2), ...
    'Jacobian', @failing_jacobian);
result = networkexp.evaluate_periodic_orbit(@vdp_rhs, [2; 0], 1, 'ode45', 100, solverOptions);
assert(result.success, 'Non-stiff solvers should ignore Jacobian metadata instead of failing.');
assert(result.hasOrbit && result.orbitCode == 2, ...
    'Successful network periodic-orbit evaluation should expose strict orbit metadata.');
end

function jacobian = finite_difference(rhs, state)
f0 = rhs(state);
numStates = numel(state);
jacobian = zeros(numel(f0), numStates);
stepSize = 1e-7;
for idx = 1:numStates
    perturbation = zeros(numStates, 1);
    perturbation(idx) = stepSize;
    jacobian(:, idx) = (rhs(state + perturbation) - f0) / stepSize;
end
end

function dydt = fhn_rhs(x, params)
theta = params{1};
gamma = params{2};
epsilon = params{3};
tau = params{4};
current = params{5};
matAdj = params{6};

numNodes = size(matAdj, 1);
dydt = zeros(2 * numNodes, 1);
for i = 1:numNodes
    dydt(i) = x(i) * (x(i) - theta) * (1 - x(i)) - x(numNodes + i) + current;
    for j = 1:numNodes
        dydt(i) = dydt(i) + matAdj(i, j) / (1 + exp(-tau * x(j)));
    end
    dydt(numNodes + i) = epsilon * (x(i) - gamma * x(numNodes + i));
end
end

function dydt = grn_rhs(y, params)
d_p = params{1};
d_E = params{2};
k_m = params{3};
epsilon = params{4};
S = params{5};
K_d = params{6};
n = params{7};
d_m = params{8};
matAdj = params{9};

numNodes = size(matAdj, 1);
dydt = zeros(2 * numNodes, 1);
for i = 1:numNodes
    p_i = y(i);
    m_i = y(i + numNodes);
    dydt(i) = m_i - d_p * p_i - d_E * p_i / (k_m + p_i + p_i^2);
    for j = 1:numNodes
        p_j = y(j);
        dydt(i) = dydt(i) + matAdj(i, j) * m_i * p_j / (1 + p_j);
    end
    dydt(numNodes + i) = epsilon * (S * K_d^n / (K_d^n + p_i^n) - d_m * m_i);
end
end

function dydt = vdp_rhs(~, y, mu)
y = y(:);
dydt = [y(2); mu * (1 - y(1)^2) * y(2) - y(1)];
end

function jacobian = failing_jacobian(~, ~)
error('test_network_solver_options:UnexpectedJacobianCall', ...
    'Jacobian should have been filtered out for non-stiff solvers.');
jacobian = [];
end
