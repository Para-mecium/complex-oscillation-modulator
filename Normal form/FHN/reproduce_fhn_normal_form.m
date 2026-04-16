function result = reproduce_fhn_normal_form(opts)
%REPRODUCE_FHN_NORMAL_FORM Build FHN normal-form data and optionally plot it.

ensure_paths();

if nargin < 1
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);
data = build_fhn_normal_form_data(cfg);

fig = [];
if cfg.make_plots
    fig = plot_fhn_normal_form(data, cfg);
end

result = struct('data', data, 'figure', fig);
end

function opts = default_options()
opts = struct();
opts.make_plots = true;
opts.show_constraint_plot = false;
end
