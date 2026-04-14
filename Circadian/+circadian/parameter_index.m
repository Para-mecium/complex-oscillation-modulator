function idx = parameter_index(pathOrModel, name)
idx = find(strcmp(pathOrModel.paramNames, name), 1, 'first');
end
