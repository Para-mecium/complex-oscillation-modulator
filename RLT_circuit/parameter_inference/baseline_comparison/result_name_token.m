function token = result_name_token(value)
token = lower(char(string(value)));
token = regexprep(token, '[^a-z0-9]+', '_');
token = regexprep(token, '^_+|_+$', '');
if isempty(token)
    token = 'unnamed';
end
end
