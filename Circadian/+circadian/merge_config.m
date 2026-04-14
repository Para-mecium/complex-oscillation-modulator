function out = merge_config(base, override)
out = base;
if nargin < 2 || isempty(override)
    return
end

fields = fieldnames(override);
for i = 1:numel(fields)
    name = fields{i};
    if isstruct(override.(name)) && isfield(out, name) && isstruct(out.(name))
        out.(name) = circadian.merge_config(out.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end
end
