function tag = format_sigma_tag(sigma)
sigmaValue = double(sigma);
if ~(isscalar(sigmaValue) && isfinite(sigmaValue))
    error('circadian:InvalidSigmaTag', ...
        'Sigma tag formatting expects a finite scalar sigma.');
end

raw = sprintf('%.12g', sigmaValue);
tag = strrep(raw, '.', 'p');
tag = strrep(tag, '+', '');
tag = strrep(tag, '-', 'm');
end
