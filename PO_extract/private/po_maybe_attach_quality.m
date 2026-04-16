function result = po_maybe_attach_quality(result, rhsFun, enabled)
if nargin < 3 || ~logical(enabled)
    return;
end

if ~isfield(result, 'diagnostics') || ~isstruct(result.diagnostics) || isempty(result.diagnostics)
    result.diagnostics = struct();
end

result.diagnostics.quality = po_evaluate_orbit_quality(result, rhsFun);
end
