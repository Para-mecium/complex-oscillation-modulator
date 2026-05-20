function perturbation = sample_perturbation(settings, nPerturbedEdges, numNodes, seed)
%SAMPLE_PERTURBATION Draw a reproducible perturbation mask and edge list.

sourceSeed = networkexp.derive_seed(seed, 11);
targetSeed = networkexp.derive_seed(seed, 29);

sourceNodes = pick_nodes( ...
    numNodes, nPerturbedEdges, settings.perturbedNodeSelectionMode, sourceSeed);
targetNodes = pick_nodes( ...
    numNodes, nPerturbedEdges, settings.inputNodeSelectionMode, targetSeed);

linearIndices = sub2ind([numNodes, numNodes], sourceNodes(:), targetNodes(:));
mask = sparse(numNodes, numNodes);
mask(linearIndices) = 1;

perturbation = struct( ...
    'seed', seed, ...
    'sourceSeed', sourceSeed, ...
    'targetSeed', targetSeed, ...
    'nPerturbedEdges', nPerturbedEdges, ...
    'sourceNodes', reshape(sourceNodes, 1, []), ...
    'targetNodes', reshape(targetNodes, 1, []), ...
    'edges', [sourceNodes(:), targetNodes(:)], ...
    'linearIndices', reshape(linearIndices, 1, []), ...
    'appliedPerturbationCount', nnz(mask), ...
    'mask', mask);
end

function nodes = pick_nodes(numNodes, numRequested, modeName, seed)
modeName = lower(char(string(modeName)));

switch modeName
    case 'prefix'
        if numRequested > numNodes
            error('networkexp:TooManyPrefixNodes', ...
                'Requested %d prefix nodes, but only %d nodes are available.', ...
                numRequested, numNodes);
        end
        nodes = 1:numRequested;

    case {'randomwithoutreplacement', 'uniformrandomwithoutreplacement'}
        if numRequested > numNodes
            error('networkexp:TooManyUniqueNodes', ...
                'Cannot sample %d unique nodes from %d nodes.', ...
                numRequested, numNodes);
        end
        previousRng = rng;
        cleanupObj = onCleanup(@() rng(previousRng));
        rng(seed, 'twister');
        nodes = randperm(numNodes, numRequested);

    case {'randomwithreplacement', 'uniformrandomwithreplacement'}
        previousRng = rng;
        cleanupObj = onCleanup(@() rng(previousRng));
        rng(seed, 'twister');
        nodes = randi(numNodes, 1, numRequested);

    otherwise
        error('networkexp:UnknownSelectionMode', ...
            'Unsupported node selection mode "%s".', modeName);
end
end
