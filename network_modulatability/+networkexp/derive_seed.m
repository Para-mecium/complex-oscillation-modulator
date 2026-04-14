function seed = derive_seed(masterSeed, varargin)
%DERIVE_SEED Deterministically derive a positive integer RNG seed.

modulus = uint64(4294967291);
state = uint64(mod(double(masterSeed), double(modulus)));
if state == 0
    state = uint64(5489);
end

for i = 1:numel(varargin)
    value = uint64(mod(double(varargin{i}), double(modulus)));
    state = mod(state * uint64(1664525) + value + uint64(1013904223 + i), modulus);
    if state == 0
        state = uint64(i);
    end
end

seed = double(state);
if seed <= 0
    seed = 1;
end
end
