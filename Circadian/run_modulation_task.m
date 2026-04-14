function result = run_modulation_task(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
result = circadian.run_fmam_task(cfg);
end
