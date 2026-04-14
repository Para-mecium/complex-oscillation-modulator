function out = matcont_base_system(varargin)
global matcont_base_system_constants

if nargin >= 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'set_constants')
    matcont_base_system_constants = varargin{2};
    out = [];
    return
end

if isempty(matcont_base_system_constants)
    cfg = default_config();
    matcont_base_system_constants = cfg.base.constants;
end

out = { ...
    @init, ...
    @fun_eval, ...
    @jacobian, ...
    @jacobianp, ...
    @hessians, ...
    @hessiansp, ...
    [], ...
    [], ...
    []};

    function dydt = fun_eval(~, state, I1, ET)
        c = current_constants();
        x = state(1);
        y = state(2);
        dydt = [ ...
            protein_production(y) * regulated_drive(I1) - c.kdx * x; ...
            c.ksy * x - c.kdy * y - c.k2 * ET * degradation_term(y)];
    end

    function [tspan, y0, options] = init
        handles = feval(@matcont_base_system);
        y0 = [0; 0];
        options = odeset( ...
            'Jacobian', handles{3}, ...
            'JacobianP', handles{4}, ...
            'Hessians', handles{5}, ...
            'HessiansP', handles{6});
        tspan = [0 10];
    end

    function jac = jacobian(~, state, I1, ET)
        c = current_constants();
        y = state(2);
        jac = [ ...
            -c.kdx, protein_production_dy(y) * regulated_drive(I1); ...
            c.ksy, -c.kdy - c.k2 * ET * degradation_term_dy(y)];
    end

    function jacp = jacobianp(~, state, I1, ET)
        c = current_constants();
        y = state(2);
        jacp = [ ...
            protein_production(y) * regulated_drive_di1(I1), 0; ...
            0, -c.k2 * degradation_term(y)];
    end

    function hess = hessians(~, state, I1, ET)
        c = current_constants();
        y = state(2);
        hess = zeros(2, 2, 2);
        hess(:, :, 1) = zeros(2, 2);
        hess(:, :, 2) = [ ...
            0, protein_production_dyy(y) * regulated_drive(I1); ...
            0, -c.k2 * ET * degradation_term_dyy(y)];
    end

    function hessp = hessiansp(~, state, I1, ET)
        c = current_constants();
        y = state(2);
        hessp = zeros(2, 2, 2);
        hessp(:, :, 1) = [ ...
            0, protein_production_dy(y) * regulated_drive_di1(I1); ...
            0, 0];
        hessp(:, :, 2) = [ ...
            0, 0; ...
            0, -c.k2 * degradation_term_dy(y)];
    end
end

function value = protein_production(y)
c = current_constants();
value = c.k1 * c.S * c.Kd^c.p / (c.Kd^c.p + y^c.p);
end

function value = protein_production_dy(y)
c = current_constants();
value = -c.k1 * c.S * c.Kd^c.p * c.p * y^(c.p - 1) / (c.Kd^c.p + y^c.p)^2;
end

function value = protein_production_dyy(y)
c = current_constants();
denominator = c.Kd^c.p + y^c.p;
value = -c.k1 * c.S * c.Kd^c.p * c.p * (c.p - 1) * y^(c.p - 2) / denominator^2 + ...
    2 * c.k1 * c.S * c.Kd^c.p * c.p^2 * y^(2 * c.p - 2) / denominator^3;
end

function value = degradation_term(y)
c = current_constants();
value = y / degradation_denominator(y, c);
end

function value = degradation_term_dy(y)
c = current_constants();
denominator = degradation_denominator(y, c);
value = (c.Km - c.KI * y^2) / denominator^2;
end

function value = degradation_term_dyy(y)
c = current_constants();
denominator = degradation_denominator(y, c);
numerator = c.Km - c.KI * y^2;
denominatorDy = 1 + 2 * c.KI * y;
value = (-2 * c.KI * y * denominator - 2 * numerator * denominatorDy) / denominator^3;
end

function value = degradation_denominator(y, c)
value = c.Km + y + c.KI * y^2;
end

function value = regulated_drive(I1)
c = current_constants();
hU = c.U * I1^2 / (c.K1 + I1^2);
value = c.bU * hU^2 / (c.KU + hU^2);
end

function value = regulated_drive_di1(I1)
c = current_constants();
hU = c.U * I1^2 / (c.K1 + I1^2);
dhU = 2 * c.U * c.K1 * I1 / (c.K1 + I1^2)^2;
value = 2 * c.bU * c.KU * hU * dhU / (c.KU + hU^2)^2;
end

function c = current_constants()
global matcont_base_system_constants
if isempty(matcont_base_system_constants)
    cfg = default_config();
    matcont_base_system_constants = cfg.base.constants;
end
c = matcont_base_system_constants;
end
