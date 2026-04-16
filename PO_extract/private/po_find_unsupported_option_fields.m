function unsupportedFields = po_find_unsupported_option_fields(opts)
%PO_FIND_UNSUPPORTED_OPTION_FIELDS Return unsupported option names from an opts struct.

if isempty(opts)
    unsupportedFields = strings(0, 1);
    return;
end

optionNames = string(fieldnames(opts));
supportedFields = po_supported_option_fields();
unsupportedFields = sort(optionNames(~ismember(optionNames, supportedFields)));
unsupportedFields = reshape(string(unsupportedFields), [], 1);
end
