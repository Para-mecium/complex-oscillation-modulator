function solverOptions = build_fhn_solver_options(params)
%BUILD_FHN_SOLVER_OPTIONS Assemble sparse solver metadata for FHN networks.

theta = params{1};
gamma = params{2};
epsilon = params{3};
tau = params{4};
matAdj = sparse(params{6});

numNodes = size(matAdj, 1);
identity = speye(numNodes);
jPattern = [spones(matAdj) + identity, identity; identity, identity];

solverOptions = struct();
solverOptions.JPattern = jPattern;
solverOptions.Jacobian = @(t, x) evaluate_fhn_jacobian(t, x, theta, gamma, epsilon, tau, matAdj); %#ok<INUSD>
end

function jacobian = evaluate_fhn_jacobian(~, x, theta, gamma, epsilon, tau, matAdj)
x = x(:);
numNodes = size(matAdj, 1);
u = x(1:numNodes);

sigma = 1 ./ (1 + exp(-tau * u));
sigmaPrime = tau * sigma .* (1 - sigma);
localDerivative = -3 * u.^2 + 2 * (1 + theta) * u - theta;

jxx = spdiags(localDerivative, 0, numNodes, numNodes) + matAdj * spdiags(sigmaPrime, 0, numNodes, numNodes);
identity = speye(numNodes);
jacobian = [jxx, -identity; epsilon * identity, -epsilon * gamma * identity];
end
