function fig = plot_fhn_normal_form(data, opts)
%PLOT_FHN_NORMAL_FORM Plot precomputed FHN normal-form curves.

ensure_paths();

if nargin < 2
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);

if cfg.show_constraint_plot
    fig = figure('Name', cfg.figure_name);
    layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    fig = figure('Name', cfg.figure_name);
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
end

am = data.curves.am;
fm = data.curves.fm;
checks = data.checks;
epsilon = data.meta.parameters.epsilon;

ax = nexttile(layout);
plot(ax, am.f11, am.f12, 'LineWidth', 1.6);
grid(ax, 'on');
xlabel(ax, 'f_{11}');
ylabel(ax, 'f_{12}');
title(ax, sprintf('AM: max |\\Delta\\omega^2| = %.3e', checks.am_omega_max_abs));

ax = nexttile(layout);
plot(ax, fm.f11, fm.f12, 'LineWidth', 1.6);
grid(ax, 'on');
xlabel(ax, 'f_{11}');
ylabel(ax, 'f_{12}');
title(ax, sprintf('FM: max |\\Delta\\chi| = %.3e', checks.fm_chi_max_abs));

if cfg.show_constraint_plot
    ax = nexttile(layout);
    theoryF12 = [min([am.f12; fm.f12]); max([am.f12; fm.f12])];
    theoryF21 = -epsilon .* theoryF12;
    plot(ax, am.f12, am.f21, 'LineWidth', 1.2);
    hold(ax, 'on');
    plot(ax, fm.f12, fm.f21, 'LineWidth', 1.2);
    plot(ax, theoryF12, theoryF21, 'k--', 'LineWidth', 1.0);
    grid(ax, 'on');
    xlabel(ax, 'f_{12}');
    ylabel(ax, 'f_{21}');
    title(ax, sprintf('Constraint: max residual = %.3e', checks.constraint_max_abs));
    legend(ax, {'AM', 'FM', 'f_{21}=-\epsilon f_{12}'}, 'Location', 'best');
end
end

function opts = default_options()
opts = struct();
opts.figure_name = 'FHN normal-form curves';
opts.show_constraint_plot = false;
end
