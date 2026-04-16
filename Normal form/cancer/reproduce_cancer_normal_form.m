function result = reproduce_cancer_normal_form(opts)
%REPRODUCE_CANCER_NORMAL_FORM Build cancer normal-form data and optionally plot it.

ensure_paths();

if nargin < 1
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);
data = build_cancer_normal_form_data(cfg);

fig = [];
if cfg.make_plots
    fig = plot_cancer_normal_form(data, cfg);
end

result = struct('data', data, 'figure', fig);
end

function opts = default_options()
opts = struct();
opts.make_plots = true;
opts.show_constraint_plot = false;
end
