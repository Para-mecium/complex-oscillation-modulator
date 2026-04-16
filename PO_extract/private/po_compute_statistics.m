function [maxVariable, minVariable, amplitude] = po_compute_statistics(orbitY)
if isempty(orbitY)
    maxVariable = [];
    minVariable = [];
    amplitude = [];
    return;
end

maxVariable = max(orbitY, [], 1);
minVariable = min(orbitY, [], 1);
amplitude = (maxVariable - minVariable) / 2;
end
