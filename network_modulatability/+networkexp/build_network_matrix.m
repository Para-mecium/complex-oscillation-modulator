function [omega, metadata] = build_network_matrix(networkSettings, netName, seed)
%BUILD_NETWORK_MATRIX Generate one weighted network with a fixed RNG seed.

packageDir = fileparts(mfilename('fullpath'));
networkDir = fileparts(packageDir);
if exist(networkDir, 'dir') == 7 && ~contains(path, networkDir)
    addpath(networkDir);
end

previousRng = rng;
cleanupObj = onCleanup(@() rng(previousRng));
rng(seed, 'twister');

N = networkSettings.N;
degAvg = networkSettings.degAvg;
K = networkSettings.K;
netName = upper(string(netName));

switch netName
    case "ER"
        adjacency = double(rand(N, N) <= degAvg / N);

    case "BA"
        ensure_generator_exists('BAgraph_dir');
        adjacency = double(BAgraph_dir(N, degAvg, networkSettings.baSeedDegree));

    case "SW"
        ensure_generator_exists('buildsw');
        adjacency = double(buildsw(N, networkSettings.swNeighborCount, networkSettings.swRewireProb));

    otherwise
        error('networkexp:UnknownNetwork', 'Unsupported network "%s".', netName);
end

adjacency = sparse(adjacency);
[rowIdx, colIdx] = find(adjacency);
edgeWeights = 0.2 + 0.8 * rand(numel(rowIdx), 1);
omega = sparse(rowIdx, colIdx, (K / N) * edgeWeights, N, N);

metadata = struct( ...
    'netName', char(netName), ...
    'seed', seed, ...
    'numNodes', N, ...
    'degreeAverage', degAvg, ...
    'couplingScale', K, ...
    'numEdges', nnz(adjacency), ...
    'density', nnz(adjacency) / numel(adjacency));
end

function ensure_generator_exists(functionName)
if exist(functionName, 'file') ~= 2
    error('networkexp:MissingGenerator', ...
        'Required generator "%s" is not on the MATLAB path.', functionName);
end
end
