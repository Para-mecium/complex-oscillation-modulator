function opts = merge_options(defaults, overrides)
%MERGE_OPTIONS Merge override fields into a defaults struct.

opts = defaults;

if nargin < 2 || isempty(overrides)
    return
end

fields = fieldnames(overrides);
for idx = 1:numel(fields)
    opts.(fields{idx}) = overrides.(fields{idx});
end
end
