function [itemsPerturb, itemsControlled, maxStep, accuracy] = build_task_inputs(taskSpec, model, solverView, derivedView, cfg)
goalOrder = taskSpec.goalOrder;
itemsPerturb = repmat(struct('prop', '', 'idx', 1, 'target', NaN), 1, numel(goalOrder));
itemsControlled = zeros(1, numel(goalOrder));
delta = zeros(1, numel(goalOrder));
accuracy = zeros(1, numel(goalOrder));

for i = 1:numel(goalOrder)
    name = goalOrder{i};
    key = lower(name);
    switch key
        case 'period'
            itemsPerturb(i).prop = 'p_Psi';
            itemsPerturb(i).idx = 1;
            itemsPerturb(i).target = taskSpec.goals.period / (2 * pi);
            delta(i) = itemsPerturb(i).target - solverView.p_Psi(1);
            accuracy(i) = target_accuracy(cfg, 'period');

        case {'amplitude', 'obsamp'}
            itemsPerturb(i).prop = 'obsAmp';
            itemsPerturb(i).idx = 1;
            itemsPerturb(i).target = goal_value(taskSpec.goals, name, 'amplitude', 'obsAmp');
            delta(i) = itemsPerturb(i).target - derivedView.obsAmp(1);
            accuracy(i) = target_accuracy(cfg, 'obsAmp', 'amplitude');

        case {'maximum', 'obsmax'}
            itemsPerturb(i).prop = 'obsMax';
            itemsPerturb(i).idx = 1;
            itemsPerturb(i).target = goal_value(taskSpec.goals, name, 'maximum', 'obsMax');
            delta(i) = itemsPerturb(i).target - derivedView.obsMax(1);
            accuracy(i) = target_accuracy(cfg, 'obsMax', 'maximum');

        case {'minimum', 'obsmin'}
            itemsPerturb(i).prop = 'obsMin';
            itemsPerturb(i).idx = 1;
            itemsPerturb(i).target = goal_value(taskSpec.goals, name, 'minimum', 'obsMin');
            delta(i) = itemsPerturb(i).target - derivedView.obsMin(1);
            accuracy(i) = target_accuracy(cfg, 'obsMin', 'minimum');

        otherwise
            idx = circadian.parameter_index(model, name);
            if isempty(idx)
                error('circadian:UnknownGoal', ...
                    'Unsupported goal ''%s''. Expected period, amplitude/obsAmp, maximum/obsMax, minimum/obsMin, or a parameter name.', ...
                    name)
            end
            itemsPerturb(i).prop = 'params';
            itemsPerturb(i).idx = idx;
            itemsPerturb(i).target = goal_value(taskSpec.goals, name);
            delta(i) = itemsPerturb(i).target - solverView.params(idx);
            accuracy(i) = target_accuracy(cfg, name, 'params');
    end

    controlName = taskSpec.controlledParams{i};
    itemsControlled(i) = circadian.parameter_index(model, controlName);
    if isempty(itemsControlled(i))
        error('circadian:UnknownControlledParameter', ...
            'Unknown controlled parameter ''%s''.', controlName)
    end
end

if isfield(taskSpec, 'maxStepScale') && ~isempty(taskSpec.maxStepScale)
    maxStepScale = double(taskSpec.maxStepScale);
else
    maxStepScale = 1;
end
if isscalar(maxStepScale)
    maxStepScale = repmat(maxStepScale, 1, numel(goalOrder));
else
    maxStepScale = reshape(maxStepScale, 1, []);
end
if numel(maxStepScale) ~= numel(goalOrder)
    error('circadian:InvalidMaxStepScale', ...
        'maxStepScale must be a scalar or one value per goal.')
end

baseStep = max(abs(delta) / cfg.fmam.targetSteps, cfg.fmam.minGoalStep);
maxStep = maxStepScale .* baseStep;
end

function value = goal_value(goals, requestedName, varargin)
aliases = [{requestedName}, varargin];
for i = 1:numel(aliases)
    alias = aliases{i};
    if isfield(goals, alias)
        value = goals.(alias);
        return
    end
end

error('circadian:MissingGoalValue', ...
    'Missing goal value for ''%s''.', requestedName)
end

function value = target_accuracy(cfg, varargin)
if isfield(cfg, 'fmam') && isfield(cfg.fmam, 'targetAccuracy')
    accuracyCfg = cfg.fmam.targetAccuracy;
    for i = 1:numel(varargin)
        name = varargin{i};
        if isfield(accuracyCfg, name)
            value = accuracyCfg.(name);
            return
        end
    end
end

value = cfg.fmam.errBound;
end
