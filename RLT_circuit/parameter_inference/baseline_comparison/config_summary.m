function summary = config_summary(config)
summary = rmfield(config, intersect(fieldnames(config), ...
    {'targetOrbit', 'targetFeatures'}));
end
