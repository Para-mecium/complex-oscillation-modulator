function tag = cache_tag(value)
if isnumeric(value)
    if ~isscalar(value)
        error('circadian:CacheTagScalarExpected', 'Numeric cache tags must be scalar.');
    end
    raw = sprintf('%.12g', value);
elseif isstring(value) || ischar(value)
    raw = char(value);
else
    error('circadian:CacheTagUnsupportedType', 'Unsupported cache tag type.');
end

raw = strrep(raw, '+', '');
raw = strrep(raw, '-', 'm');
raw = strrep(raw, '.', 'p');
raw = regexprep(raw, '[^A-Za-z0-9_]+', '_');
raw = regexprep(raw, '_+', '_');
raw = regexprep(raw, '^_+|_+$', '');
if isempty(raw)
    raw = 'value';
end

tag = raw;
end
