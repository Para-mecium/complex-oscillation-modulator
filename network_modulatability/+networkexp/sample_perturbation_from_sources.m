function perturbation = sample_perturbation_from_sources(settings, sourceNodes, numNodes, seed)
%SAMPLE_PERTURBATION_FROM_SOURCES Draw targets for a fixed source-node group.

targetSeed = networkexp.derive_seed(seed, 29);
targetNodes = pick_target_nodes( ...
    numNodes, numel(sourceNodes), settings.inputNodeSelectionMode, targetSeed);

linearIndices = sub2ind([numNodes, numNodes], sourceNodes(:), targetNodes(:));
mask = sparse(numNodes, numNodes);
mask(linearIndices) = 1;

perturbation = struct( ...
    'seed', seed, ...
    'sourceSeed', [], ...
    'targetSeed', targetSeed, ...
    'nPerturbedEdges', numel(sourceNodes), ...
    'sourceNodes', reshape(sourceNodes, 1, []), ...
    'targetNodes', reshape(targetNodes, 1, []), ...
    'edges', [sourceNodes(:), targetNodes(:)], ...
    'linearIndices', reshape(linearIndices, 1, []), ...
    'appliedPerturbationCount', nnz(mask), ...
    'mask', mask);
end

function nodes = pick_target_nodes(numNodes, numRequested, modeName, seed)
modeName = lower(char(string(modeName)));

switch modeName
    case 'prefix'
        nodes = 1:numRequested;

    case {'randomwithoutreplacement', 'uniformrandomwithoutreplacement'}
        previousRng = rng;
        cleanupObj = onCleanup(@() rng(previousRng));
        rng(seed, 'twister');
        nodes = randperm(numNodes, numRequested);

    otherwise
        error('networkexp:UnknownSelectionMode', ...
            'Unsupported target-node selection mode "%s".', modeName);
end
end
