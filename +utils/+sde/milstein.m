function [lt, X_M] = milstein(FuncDrift, FuncDiffusion, diff_FuncDiffusion, T, N, x0, sigma, noiseClass)
if ~iscell(FuncDiffusion)
    FuncDiffusion = {FuncDiffusion};
end
if ~iscell(diff_FuncDiffusion)
    diff_FuncDiffusion = {diff_FuncDiffusion};
end

x0 = reshape(double(x0), [], 1);
noiseNum = numel(FuncDiffusion);
sigma = normalize_sigma(sigma, noiseNum);
noiseClass = normalize_noise_class(noiseClass, noiseNum);

lt = linspace(0, T, N + 1);
dt = T / N;
X_M = zeros(numel(x0), N + 1);
X_M(:, 1) = x0;
ouState = zeros(noiseNum, N + 1);

for j = 1:N
    xj = X_M(:, j);
    xNext = xj + FuncDrift(lt(j), xj) * dt;

    for i = 1:noiseNum
        Fi = FuncDiffusion{i};
        switch lower(noiseClass{i})
            case 'w'
                diffFi = diff_FuncDiffusion{i};
                dW = sigma(i) * sqrt(dt) * randn;
                FiEval = Fi(lt(j), xj);
                xNext = xNext + FiEval * dW;
                xNext = xNext + 0.5 * diffFi(lt(j), xj) * FiEval * (dW^2 - sigma(i)^2 * dt);
            case 'o'
                % Match the legacy OU-driven Euler state update used by the paper code.
                xNext = xNext + Fi(lt(j), xj) * ouState(i, j) * dt;
            otherwise
                error('utils:sde:milstein:UnknownNoiseClass', ...
                    'Unknown NoiseClass "%s".', noiseClass{i});
        end
    end

    X_M(:, j + 1) = xNext;

    for i = 1:noiseNum
        if strcmpi(noiseClass{i}, 'o')
            deltaW = sqrt(dt) * randn;
            ouState(i, j + 1) = ouState(i, j) ...
                - ouState(i, j) * dt ...
                + sigma(i) * deltaW ...
                + 0.5 * sigma(i)^2 * (deltaW^2 - dt);
        end
    end
end
end

function sigmaVec = normalize_sigma(sigma, expectedNum)
sigmaVec = reshape(double(sigma), 1, []);
if isscalar(sigmaVec)
    sigmaVec = repmat(sigmaVec, 1, expectedNum);
elseif numel(sigmaVec) ~= expectedNum
    error('utils:sde:milstein:InvalidSigma', ...
        'Expected sigma to be scalar or length %d.', expectedNum);
end
end

function noiseClassOut = normalize_noise_class(noiseClass, expectedNum)
if ischar(noiseClass) || isstring(noiseClass)
    noiseClassOut = cellstr(noiseClass);
else
    noiseClassOut = noiseClass;
end

if isscalar(noiseClassOut)
    noiseClassOut = repmat(noiseClassOut, 1, expectedNum);
end
if numel(noiseClassOut) ~= expectedNum
    error('utils:sde:milstein:InvalidNoiseClass', ...
        'Expected noiseClass to be scalar or length %d.', expectedNum);
end

noiseClassOut = cellfun(@char, noiseClassOut, 'UniformOutput', false);
end
