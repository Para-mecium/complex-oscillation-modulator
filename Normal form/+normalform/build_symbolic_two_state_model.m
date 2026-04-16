function model = build_symbolic_two_state_model(spec)
%BUILD_SYMBOLIC_TWO_STATE_MODEL Build cached A/H/T handles for a shifted two-state system.

arguments
    spec struct
end

if ~isfield(spec, 'rhs_builder') || ~isa(spec.rhs_builder, 'function_handle')
    error('normalform:InvalidSpec', ...
        'spec.rhs_builder must be a function handle.');
end

require_symbolic_toolbox();

cacheKey = '';
if isfield(spec, 'cache_key') && ~isempty(spec.cache_key)
    cacheKey = string(spec.cache_key);
end

persistent cache
if isempty(cache)
    cache = struct('key', {{}}, 'value', {{}});
end

if strlength(cacheKey) > 0
    matchIdx = find(strcmp(cache.key, char(cacheKey)), 1, 'first');
    if ~isempty(matchIdx)
        model = cache.value{matchIdx};
        return
    end
end

syms u1 u2 f11 f12 real

rhs = spec.rhs_builder(u1, u2, f11, f12);
rhs = rhs(:);
if numel(rhs) ~= 2
    error('normalform:InvalidSpec', ...
        'spec.rhs_builder must return a two-element vector field.');
end

vars = [u1; u2];
zeroPoint = [0, 0];

linearTerm = subs(jacobian(rhs, vars), [u1, u2], zeroPoint);
quadraticTensor = sym(zeros(2, 2, 2));
cubicTensor = sym(zeros(2, 2, 2, 2));

for i = 1:2
    for j = 1:2
        for k = 1:2
            quadraticTensor(i, j, k) = subs(diff(rhs(i), vars(j), vars(k)), [u1, u2], zeroPoint);
            for l = 1:2
                cubicTensor(i, j, k, l) = subs(diff(rhs(i), vars(j), vars(k), vars(l)), [u1, u2], zeroPoint);
            end
        end
    end
end

model = struct();
model.A_function = matlabFunction(linearTerm, 'Vars', {f11, f12});
model.H_function = matlabFunction(quadraticTensor, 'Vars', {f11, f12});
model.T_function = matlabFunction(cubicTensor, 'Vars', {f11, f12});
model.cache_key = char(cacheKey);

if strlength(cacheKey) > 0
    cache.key{end + 1, 1} = char(cacheKey);
    cache.value{end + 1, 1} = model;
end
end

function require_symbolic_toolbox()
if exist('sym', 'file') ~= 2 || ~license('test', 'Symbolic_Toolbox')
    error('normalform:MissingSymbolicToolbox', ...
        'Automatic normal-form reconstruction requires Symbolic Math Toolbox.');
end
end
