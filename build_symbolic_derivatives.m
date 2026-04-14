function derivatives = build_symbolic_derivatives(system, observables, dimParams)
%BUILD_SYMBOLIC_DERIVATIVES Precompute analytic derivative handles for FMAM_ODE.
%   This helper keeps Symbolic Toolbox work outside FMAM_ODE so the solver
%   only consumes prepared derivatives during construction.

    if ~iscell(system) || isempty(system)
        error('build_symbolic_derivatives:InvalidSystem', ...
            'system must be a nonempty cell array of function handles.')
    end
    system = reshape(system, 1, []);

    if isempty(observables)
        observables = {};
    elseif iscell(observables)
        observables = reshape(observables, 1, []);
    else
        error('build_symbolic_derivatives:InvalidObservables', ...
            'observables must be empty or a cell array of function handles.')
    end

    if ~isscalar(dimParams) || dimParams < 1 || dimParams ~= floor(dimParams)
        error('build_symbolic_derivatives:InvalidParameterDimension', ...
            'dimParams must be a positive integer scalar.')
    end

    dimVar = numel(system);
    dimObs = numel(observables);

    derivatives = struct;
    derivatives.var = struct;
    derivatives.obs = [];
    derivatives.obs2 = [];

    parameters_sym = sym('parameter', [1, dimParams]);
    variable_sym = sym('variable', [1, dimVar]);

    for i = 1:dimVar
        fi = system{i};
        J = jacobian(fi(variable_sym, parameters_sym), [variable_sym, parameters_sym]);
        for j = 1:(dimVar + dimParams)
            derivatives.var(i, j).function = matlabFunction(J(1, j), ...
                'vars', {variable_sym, parameters_sym});
        end
    end

    if dimObs == 0
        return
    end

    derivatives.obs = struct;
    derivatives.obs2 = struct;
    for i = 1:dimObs
        fi = observables{i};
        J = jacobian(fi(variable_sym), variable_sym);
        for j = 1:dimVar
            firstDerivative = matlabFunction(J(1, j), 'vars', {variable_sym});
            derivatives.obs(i, j).function = firstDerivative;

            JJ = jacobian(firstDerivative(variable_sym), variable_sym);
            for k = 1:dimVar
                derivatives.obs2(i, j, k).function = matlabFunction(JJ(1, k), ...
                    'vars', {variable_sym});
            end
        end
    end
end
