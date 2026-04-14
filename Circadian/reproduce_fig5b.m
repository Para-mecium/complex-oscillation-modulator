function result = reproduce_fig5b(cfg)
if nargin < 1
    cfg = struct();
end

circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);
curveData = build_fig5b_data(cfg);
markData = build_fig5b_marks(cfg);

fig = plot_data(curveData, markData, cfg);
circadian.save_figure(fig, fullfile(cfg.io.figureDir, 'fig5b'));
data = struct();
data.curves = curveData.curves;
data.curveFiles = curveData.curveFiles;
data.markResults = markData.markResults;
data.markFiles = markData.markFiles;
data.markAmplitudes = markData.markAmplitudes;
result = struct('data', data, 'figure', fig);
end

function fig = plot_data(curveData, markData, cfg)
fig = figure('Color', 'w', 'Position', [100, 100, 640, 520]);
ax = axes(fig);
hold(ax, 'on');
colormap(ax, cfg.plot.amplitudeColormap);

ampValues = [];
for i = 1:numel(curveData.curves)
    ampValues = [ampValues; curve_amplitude_values(curveData.curves{i})]; 
end
if isempty(ampValues)
    ampLimits = [0, 1];
else
    ampLimits = [min(ampValues), max(ampValues)];
end

for i = 1:numel(curveData.curves)
    curve = curveData.curves{i};
    plot_curve_branches(ax, cfg, curve);
    labelAnchor = curve_label_anchor(curve);
    if ~isempty(labelAnchor)
        text(ax, cfg.plot.kdScale * labelAnchor.Kd + 0.02, labelAnchor.AT, ...
            sprintf('T=%.1f', curve.targetPeriod), ...
            'FontSize', 11, 'Color', [0.12, 0.32, 0.68]);
    end
end

for i = 1:numel(markData.markResults)
    mark = markData.markResults{i};
    scatter(ax, cfg.plot.kdScale * mark.params(1), mark.params(2), 120, mark.measures.amplitude, ...
        'filled', 'MarkerEdgeColor', [0.25, 0.25, 0.25], 'LineWidth', 1.0);
end

grid(ax, 'on');
clim(ax, ampLimits);
cb = colorbar(ax);
cb.Label.String = 'Amplitude of P_{tot}';
xlabel(ax, 'K_d (\times 10^{-4})');
ylabel(ax, 'A_T (a.u.)');
xlim(ax, cfg.plot.kdScale * cfg.fig5b.KdAxis);
ylim(ax, cfg.fig5b.ATAxis);
title(ax, 'Fig. 5b: Iso-period curves');
end

function plot_curve_branches(ax, cfg, curve)
components = curve_plot_components(curve);
for i = 1:numel(components)
    component = components{i};
    if isempty(component.Kd)
        continue
    end
    circadian.plot_gradient_curve(ax, ...
        cfg.plot.kdScale * component.Kd, component.AT, component.obsAmplitude, 3);
end
end

function components = curve_plot_components(curve)
lower = join_branches_for_plot(curve.branches.lowerLeft, curve.branches.lowerRight);
upper = join_branches_for_plot(curve.branches.upperLeft, curve.branches.upperRight);

if isempty(lower.Kd) && isempty(upper.Kd)
    components = {lower};
    return
end

if isempty(lower.Kd)
    components = {upper};
    return
end

if isempty(upper.Kd)
    components = {lower};
    return
end

components = {stitch_components_at_one_end(lower, upper)};
end

function joined = join_branches_for_plot(branchA, branchB)
joined = struct('Kd', zeros(0, 1), 'AT', zeros(0, 1), 'obsAmplitude', zeros(0, 1));

[KdA, ATA, ampA] = branch_plot_vectors(branchA);
[KdB, ATB, ampB] = branch_plot_vectors(branchB);

if isempty(KdA) && isempty(KdB)
    return
end

if isempty(KdA)
    joined.Kd = KdB;
    joined.AT = ATB;
    joined.obsAmplitude = ampB;
    return
end

if isempty(KdB)
    joined.Kd = KdA;
    joined.AT = ATA;
    joined.obsAmplitude = ampA;
    return
end

if abs(KdA(end) - KdB(1)) < 1e-12 && abs(ATA(end) - ATB(1)) < 1e-12
    KdB(1) = [];
    ATB(1) = [];
    ampB(1) = [];
end

joined = append_plot_vectors(KdA, ATA, ampA, KdB, ATB, ampB);
end

function [Kd, AT, amp] = branch_plot_vectors(branch)
Kd = branch.Kd(:);
AT = branch.AT(:);
amp = branch.obsAmplitude(:);
end

function component = reverse_component(component)
component.Kd = flipud(component.Kd(:));
component.AT = flipud(component.AT(:));
component.obsAmplitude = flipud(component.obsAmplitude(:));
end

function joined = append_components(componentA, componentB)
joined = append_plot_vectors( ...
    componentA.Kd, componentA.AT, componentA.obsAmplitude, ...
    componentB.Kd, componentB.AT, componentB.obsAmplitude);
end

function stitched = stitch_components_at_one_end(componentA, componentB)
candidates = { ...
    struct('A', componentA, 'B', componentB), ...
    struct('A', componentA, 'B', reverse_component(componentB)), ...
    struct('A', reverse_component(componentA), 'B', componentB), ...
    struct('A', reverse_component(componentA), 'B', reverse_component(componentB))};

bestIdx = 1;
bestGap = endpoint_gap(candidates{1}.A, candidates{1}.B);
for i = 2:numel(candidates)
    gap = endpoint_gap(candidates{i}.A, candidates{i}.B);
    if gap < bestGap
        bestGap = gap;
        bestIdx = i;
    end
end

stitched = append_components(candidates{bestIdx}.A, candidates{bestIdx}.B);
end

function joined = append_plot_vectors(KdA, ATA, ampA, KdB, ATB, ampB)
joined = struct('Kd', KdA(:), 'AT', ATA(:), 'obsAmplitude', ampA(:));

if isempty(KdA)
    joined.Kd = KdB(:);
    joined.AT = ATB(:);
    joined.obsAmplitude = ampB(:);
    return
end

if isempty(KdB)
    return
end

if abs(KdA(end) - KdB(1)) < 1e-12 && abs(ATA(end) - ATB(1)) < 1e-12
    KdB(1) = [];
    ATB(1) = [];
    ampB(1) = [];
end

joined.Kd = [joined.Kd; KdB(:)];
joined.AT = [joined.AT; ATB(:)];
joined.obsAmplitude = [joined.obsAmplitude; ampB(:)];
end

function gap = endpoint_gap(componentA, componentB)
if isempty(componentA.Kd) || isempty(componentB.Kd)
    gap = inf;
    return
end

gap = norm([componentA.Kd(end) - componentB.Kd(1), componentA.AT(end) - componentB.AT(1)], inf);
end

function values = curve_amplitude_values(curve)
values = [];
branchNames = {'lowerLeft', 'lowerRight', 'upperLeft', 'upperRight'};
for i = 1:numel(branchNames)
    branch = curve.branches.(branchNames{i});
    values = [values; branch.obsAmplitude(:)]; %#ok<AGROW>
end
end

function anchor = curve_label_anchor(curve)
anchor = [];
branchNames = {'upperRight', 'upperLeft', 'lowerRight', 'lowerLeft'};
for i = 1:numel(branchNames)
    branch = curve.branches.(branchNames{i});
    if isempty(branch.Kd)
        continue
    end
    anchor = struct('Kd', branch.Kd(end), 'AT', branch.AT(end));
    return
end
end
