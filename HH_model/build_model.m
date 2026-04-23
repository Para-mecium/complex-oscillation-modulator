function dydt = build_model(~, y, I, p, N)
%BUILD_MODEL RHS for an N-neuron coupled Hodgkin-Huxley model.
% State layout:
%   y = [V_1; ...; V_N; m_1; ...; m_N; h_1; ...; h_N; n_1; ...; n_N]

validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'N');
validateattributes(y, {'numeric'}, {'column', 'numel', 4 * N}, mfilename, 'y');
validateattributes(I, {'numeric'}, {'column', 'numel', N}, mfilename, 'I');

couplingType = "gap";
if isfield(p, 'couplingType')
    couplingType = lower(string(p.couplingType));
end

requiredFields = {'C', 'gNa', 'ENa', 'gK', 'EK', 'gL', 'EL'};
for k = 1:numel(requiredFields)
    fieldName = requiredFields{k};
    if ~isfield(p, fieldName)
        error('build_model:MissingField', 'Missing required parameter p.%s.', fieldName);
    end
end

V = y(1:N);
m = y(N + 1:2 * N);
h = y(2 * N + 1:3 * N);
n = y(3 * N + 1:4 * N);

[am, bm, ah, bh, an, bn] = get_rates(V);

if N == 1
    I_cpl = 0;
else
    switch couplingType
        case "gap"
            validate_coupling_field(p, 'G', [N, N], mfilename);
            G = p.G;
            Goff = G - diag(diag(G));
            I_cpl = Goff * V - sum(Goff, 2) .* V;
        case "synapse"
            validate_coupling_field(p, 'G', [N, N], mfilename);
            validate_coupling_field(p, 'tau', [N, 1], mfilename);
            validate_coupling_field(p, 'Vstar', [N, 1], mfilename);
            validate_coupling_field(p, 'Esyn', [N, 1], mfilename);

            I_cpl = zeros(N, 1);
            for i = 1:N
                for j = 1:N
                    if j ~= i
                        sj = 1 ./ (1 + exp(-p.tau(j) .* (V(j) - p.Vstar(j))));
                        I_cpl(i) = I_cpl(i) + p.G(i, j) .* sj .* (p.Esyn(i) - V(i));
                        % I_cpl(i) = I_cpl(i) + p.G(i, j) .* sj;
                    end
                end
            end
        otherwise
            error('build_model:InvalidCouplingType', ...
                'Unsupported p.couplingType: %s.', couplingType);
    end
end

I_ion = p.gNa .* (m .^ 3) .* h .* (V - p.ENa) ...
      + p.gK .* (n .^ 4) .* (V - p.EK) ...
      + p.gL .* (V - p.EL);

dVdt = (I - I_ion + I_cpl) ./ p.C;
dmdt = am .* (1 - m) - bm .* m;
dhdt = ah .* (1 - h) - bh .* h;
dndt = an .* (1 - n) - bn .* n;

dydt = [dVdt; dmdt; dhdt; dndt];
end

function [am, bm, ah, bh, an, bn] = get_rates(V)
am = 0.1 .* (V + 40) ./ (1 - exp(-(V + 40) ./ 10));
bm = 4 .* exp(-(V + 65) ./ 18);
ah = 0.07 .* exp(-(V + 65) ./ 20);
bh = 1 ./ (1 + exp(-(V + 35) ./ 10));
an = 0.01 .* (V + 55) ./ (1 - exp(-(V + 55) ./ 10));
bn = 0.125 .* exp(-(V + 65) ./ 80);

am(isnan(am)) = 1;
an(isnan(an)) = 0.1;
end

function validate_coupling_field(p, fieldName, fieldSize, functionName)
if ~isfield(p, fieldName)
    error('build_model:MissingField', 'Missing required parameter p.%s.', fieldName);
end
validateattributes(p.(fieldName), {'numeric'}, {'size', fieldSize}, functionName, ['p.' fieldName]);
end
