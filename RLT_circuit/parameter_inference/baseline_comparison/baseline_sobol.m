function result = baseline_sobol(config, seed)
rng(seed, 'twister');

result = empty_result('sobol', config, seed);
lowerBound = config.activeLowerBound;
upperBound = config.activeUpperBound;
span = upperBound - lowerBound;
numSamples = config.evalBudget;
activeDim = numel(lowerBound);

unitSamples = simple_sobol_points(numSamples, activeDim);
unitSamples = mod(unitSamples + rand(1, activeDim), 1);

runTimer = tic;
for evalIdx = 1:numSamples
    activeParams = lowerBound + span .* unitSamples(evalIdx, :);
    evalOut = evaluate_candidate(activeParams, config);
    elapsedTime = toc(runTimer);
    result = record_candidate(result, evalIdx, evalOut, config, elapsedTime);
    log_progress(result, config, evalIdx, elapsedTime);
end
result.runtime = toc(runTimer);
result = trim_result(result);
end

function points = simple_sobol_points(numSamples, dim)
if dim > 4
    error('baseline_sobol:UnsupportedDimension', ...
        'The lightweight Sobol implementation supports up to 4 dimensions.');
end

maxBits = 32;
directions = zeros(dim, maxBits, 'uint32');

directions(1, :) = uint32(2.^(31:-1:0));
primitive = [ ...
    struct('s', 1, 'a', 0, 'm', 1); ...
    struct('s', 2, 'a', 1, 'm', [1 3]); ...
    struct('s', 3, 'a', 1, 'm', [1 3 1])];

for dimIdx = 2:dim
    p = primitive(dimIdx - 1);
    s = p.s;
    a = p.a;
    m = p.m;

    for bitIdx = 1:s
        directions(dimIdx, bitIdx) = uint32(m(bitIdx) * 2^(32 - bitIdx));
    end

    for bitIdx = (s + 1):maxBits
        value = bitxor(directions(dimIdx, bitIdx - s), ...
            bitshift(directions(dimIdx, bitIdx - s), -s));

        for k = 1:(s - 1)
            if bitget(a, s - k)
                value = bitxor(value, directions(dimIdx, bitIdx - k));
            end
        end

        directions(dimIdx, bitIdx) = value;
    end
end

points = zeros(numSamples, dim);
x = zeros(1, dim, 'uint32');
for sampleIdx = 1:numSamples
    bitIdx = trailing_zero_count(uint32(sampleIdx)) + 1;
    for dimIdx = 1:dim
        x(dimIdx) = bitxor(x(dimIdx), directions(dimIdx, bitIdx));
        points(sampleIdx, dimIdx) = double(x(dimIdx)) / 2^32;
    end
end
end

function count = trailing_zero_count(value)
count = 0;
while bitand(value, uint32(1)) == 0
    count = count + 1;
    value = bitshift(value, -1);
end
end
