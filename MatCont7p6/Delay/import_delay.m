function varargout = import_delay(varargin)
% SYSTEM Application M-file for system.fig
%    FIG = SYSTEM launch system GUI.
%    SYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 08-Feb-2024 23:41:07
global gds oldgds path_sys MC driver_window;

if nargin == 0 ||((nargin ==1)&&(strcmp(varargin{1},'new'))) % LAUNCH GUI
    h=gcbo;
    arg=get(h,'Tag');
    fig = openfig(mfilename,'reuse', 'invisible');
 % Use system color scheme for figure:
set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

   if strcmp(arg,'import_delay')||((nargin ==1)&&strcmp(varargin{1},'new'))
        if ~isempty(MC.starter),gds.open.figuur=1;else gds.open.figuur=0;end
        delete(MC.starter);MC.starter=[];
        % continuer-window open?
        if ~isempty(MC.continuer), gds.open.continuer=1;else gds.open.continuer=0;end
        delete(MC.continuer);MC.continuer=[];
        % numeric window open?   
        if ~isempty(MC.numeric_fig), gds.open.numeric_fig=1;else gds.open.numeric_fig=0;end
        close(MC.numeric_fig);MC.numeric_fig=[];
        %2D-plot open      
        if ~isempty(MC.D2), gds.open.D2=size(MC.D2);else gds.open.D2=0;end
        close(MC.D2);MC.D2=[];
        %3D-plot open      
        if ~isempty(MC.D3), gds.open.D3=size(MC.D3);else gds.open.D3=0;end
        close(MC.D3);MC.D3=[];%   
        %PRC-plot open      
        if ~isempty(MC.PRC), gds.open.PRC=size(MC.PRC);else gds.open.PRC=0;end
        close(MC.PRC);MC.PRC=[];%    
        %dPRC-plot open      
        if ~isempty(MC.dPRC), gds.open.dPRC=size(MC.dPRC);else gds.open.dPRC=0;end
        close(MC.dPRC);MC.dPRC=[];%   
        if ~isempty(MC.integrator),gds.open.integrator=1;else gds.open.integrator=0;end;
        delete(MC.integrator);MC.integrator=[];
        if isfield(gds,'der')
            oldgds=gds;
            init;
        else
            init;
            oldgds=gds;
        end
    end
    % Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
   guidata(fig, handles);
    load_system(handles);
    if nargout > 0
		varargout{1} = fig;
    end

    % http://undocumentedmatlab.com/blog/customizing-listbox-editbox-scrollbars/
    try % might not work on older versions of matlab
        jScrollPane = findjobj(handles.sys);
        set(jScrollPane, 'HorizontalScrollBarPolicy', ...
            jScrollPane.java.HORIZONTAL_SCROLLBAR_ALWAYS);
        cbStr = sprintf('set(gcbo,''HorizontalScrollBarPolicy'',%d)', ...
            jScrollPane.java.HORIZONTAL_SCROLLBAR_ALWAYS);
        hjScrollPane = handle(jScrollPane,'CallbackProperties');
        set(hjScrollPane,'ComponentResizedCallback',cbStr);
        jViewPort = jScrollPane.getViewport;
        jEditbox = jViewPort.getComponent(0);
        jEditbox.setWrapping(false);
    end

    movegui(fig, 'center');
    fig.Position(3) = fig.Position(3)*1.4;
    set(fig, 'visible', 'on');
    gds.ok = false;
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch ME
        fprintf(getReport(ME));
        %disp(ME)
		errordlg(ME.message,mfilename);
%         delete(driver_window);  
        global waithndl;
        delete(waithndl);
        waithndl = [];
    end

end

% --------------------------------------------------------------------
function ok_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok.
global gds hds path_sys MC;
filterspaces = @(z) z(z ~= ' ');
string_jac='';string_jacp='';string_hess='';string_hessp='';string_tensor3='';string_tensor4='';string_tensor5='';
if isempty(gds.system)
    errordlg('You have to give the system a name','Name');
    return
end

original_cor=filterspaces(get(handles.coordinates,'String'));
original_pa=filterspaces(get(handles.parameters,'String'));
original_equations=deblank(get(handles.sys,'String'));
[equations, cor, pa] = renameforsym(original_equations, original_cor, original_pa);

gds.delay.equations = equations;
gds.parameters = toGdsStruct(pa);
gds.delay.coordinates = toGdsStruct(cor);
par = pa;
if isempty(par)
    button = questdlg('The system has no parameters! Do you want to continue? ',...
    'Parameters','Yes','No','No');
    if strcmp(button,'No')
        return
    end
end
if (~isempty(par))
   par=strcat(',',par);
end
t=get(handles.time,'String');
if (~isempty(t))
    t=strcat(t,',');
end

try
    [num_delays, delays_string] = extract_delays;
catch ME
    errordlg(ME.message,'Error')
    return
end
try
    check_delays(delays_string);
catch
    errordlg('Some delay expressions are not valid! The evaluation of delay expressions (expr in [t-expr] or [t+expr]) should result in a numeric value; the delays may depend on the parameters but not on the time variable or on the coordinates.','Error')
    return
end
[unit_nodes, weights, unit_DM] = cheb(gds.delay.CollocationDegree);
[unit_ccweights, unit_ccnodes] = clencurt(gds.delay.QuadratureDegree);
try
    [string_sys_temp,string_sys_eq]=replace_sys_input(num_delays);
catch ME
    errordlg(ME.message,'Error')
    return
end
[ind_d,ind_r,dim_d,dim_r,ind_d_all,ind_r_all] = create_indices;
gds.dim = dim_d+gds.delay.CollocationDegree*gds.delay.dim;
if size(string_sys_temp,1) > 0
    string_sys_temp = [string_sys_temp; {''}];
end
% unused lines will be commented
% ("_d"=differential, "_r"=renewal, "_c"=coupled)
if dim_r == 0
    % only differential equations
    c_d = '';
    c_r = '%';
    c_c = '%';
    c_d_or_c = '';
    c_r_or_c = '%';
else
    c_d = '%';
    c_r_or_c = '';
    if dim_d == 0
        % only renewal equations
        c_r = '';
        c_c = '%';
        c_d_or_c = '%';
    else
        % coupled differential and renewal equations
        c_r = '%';
        c_c = '';
        c_d_or_c = '';
    end
end
string_sys_dimensions = [{
    sprintf('M = %d; %% degree of collocation polynomials', gds.delay.CollocationDegree);
    '';
    sprintf('dim = %d; %% total dimension of the delay equation', gds.delay.dim);
    sprintf('dim_d = %d; %% dimension of the differential part', dim_d);
    sprintf('dim_r = %d; %% dimension of the renewal part', dim_r);
    [c_d_or_c 'ind_d = [' print_matrix_oneline_string(ind_d) ']; % indices of differential equations'];
    [c_r_or_c 'ind_r = [' print_matrix_oneline_string(ind_r) ']; % indices of renewal equations'];
    '% indices of all components of differential part in kmrgd:';
    [c_d_or_c 'ind_d_all = [' print_matrix_oneline_string(ind_d_all) '];'];
    '% indices of all components of renewal part in kmrgd:';
    [c_r_or_c 'ind_r_all = [' print_matrix_oneline_string(ind_r_all) '];'];
    '';
    }];
string_sys_nodes_etc = [{
    '';
    '% Chebyshev type II nodes (extrema) in [-1, 1]:';
    ['unit_nodes = [' print_matrix_oneline_string(unit_nodes) '];'];
    'nodes = (unit_nodes-1)*tau_max/2; % Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]';
    '% weights for barycentric interpolation on Chebyshev type II nodes (extrema) (scaling-independent):';
    ['weights = [' print_matrix_oneline_string(weights) '];'];
    '% differentiation matrix for Chebyshev type II nodes (extrema) in [-1, 1]:';
    'unit_DM = [';
    };
    print_matrix_multiline_cell(unit_DM);
    {
    '];';
    'DM = unit_DM*2/tau_max; % differentiation matrix for Chebyshev type II nodes (extrema) in the delay interval [-tau_max, 0]';
    '';
    }];
string_sys_kmrgd_r_der = [{
    '% values of the derivative of the renewal part of the (integrated) state at all the nodes';
    '%     (M+1 blocks of length dim_r):';
    [c_r 'kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd,dim_r,M))),(M+1)*dim_r,1); % only renewal equations'];
    [c_c 'kmrgd_r_der = reshape(transpose(DM(:,2:end)*transpose(reshape(kmrgd(ind_r_all),dim_r,M))),(M+1)*dim_r,1); % coupled differential and renewal equations'];
    '% equivalent to kron(DM(:,2:end),eye(dim_r))*kmrgd(ind_r_all)';
    '';
    }];
string_sys_dydt = [
    string_sys_dimensions;
    {
    ['delays = ' delays_string ';'];
    'if any(delays < 0)';
    '    % the equation is advanced instead of delayed: we need to stop';
    '    dydt = NaN(size(kmrgd));';
    '    return';
    'else';
    '    tau_max = max(delays);';
    '    if tau_max == 0';
    '        % the equation is not delayed: we need to stop';
    '        dydt = NaN(size(kmrgd));';
    '        return';
    '    end';
    'end';
    };
    string_sys_nodes_etc;
    {
    '% kmrgd has length dim_d+M*dim';
    '% kmrgd(i), i in 1:dim_d = value of the i-th differential component of the state';
    '%     (i.e. the ind_d(i)-th component of the full state) at the current time';
    '% kmrgd(dim_d+(m-1)*dim+i), m in 1:M, i in 1:dim_d = value of the i-th component';
    '%     of the (full) state at the m-th node in the past';
    '% kmrgd([i,(dim_d+ind_d(i)):dim:(dim_d+M*dim))]), i in 1:dim_d = values of the';
    '%     i-th differential component of the state (i.e. the ind_d(i)-th component of';
    '%     the full state) at all nodes (length M+1); equivalent to kmrgd(ind_d_all(i:dim_d:end))';
    '% kmrgd(dim_d+(m-1)*dim+(1:dim)), m in 1:M = values of the (full) state at the m-th node in the past';
    '% kmrgd(ind_d_all) = values of the differential components of the state at all the nodes';
    '%     (M+1 blocks of length dim_d)';
    '% kmrgd(ind_r_all) = values of the (integrated) renewal components of the state at';
    '%     the nodes in the past (value at current time is always 0) (M blocks of length dim_r)';
    '';
    };
    string_sys_kmrgd_r_der;
    {
    '% kmrgd_r_der(i), i in 1:dim_r = value of the derivative of the i-th (integrated)';
    '%     renewal component of the state (i.e. the ind_r(i)-th component of the full state)';
    '%     at the current time (i.e. true value of the i-th component)';
    '% kmrgd_r_der(m*dim_r+i), m in 1:M, i in 1:dim_r = value of the same at m-th node in the past';
    '% kmrgd_r_der(i:dim_r:(M+1)*dim_r), i in 1:dim_r = values of the same at all nodes (length M+1)';
    '';
    };
    string_sys_temp;
    {
    '% renewal part of the right-hand side:';
    [c_r_or_c 'F_r = ['];
    };
    string_sys_eq(ind_r,1);
    {
    [c_r_or_c '];'];
    '% F_r extended with 0 for the differential part:';
    '% (ind_r_all-dim_d removes the offset in the indices due to the first block in kmrgd)';
    [c_r 'F_r0 = repmat(F_r,M,1); % only renewal equations'];
    [c_c 'F_r0 = zeros(M*dim,1); F_r0(ind_r_all-dim_d) = repmat(F_r,M,1); % coupled differential and renewal equations'];
    '';
    '% kmrgd with added zeros for the renewal components at current time (M+1 blocks of length dim):';
    [c_d 'kmrgd0 = kmrgd; % only differential equations'];
    [c_r 'kmrgd0 = [zeros(dim_r,1); kmrgd]; % only renewal equations'];
    [c_c 'kmrgd0 = zeros((M+1)*dim,1); kmrgd0([ind_d,dim+1:end]) = kmrgd; % coupled differential and renewal equations'];
    '';
    'dydt=['
    };
    string_sys_eq(ind_d,1);
    {
    [c_d 'reshape(transpose(DM(2:end,:)*transpose(reshape(kmrgd0,dim,M+1))),M*dim,1); % only differential equations'];
    [c_r_or_c 'reshape(transpose(DM(2:end,:)*transpose(reshape(kmrgd0,dim,M+1))),M*dim,1)-F_r0; % only/also renewal equations'];
    '];';
    '% in the last line, reshape(...) is equivalent to kron(DM(2:end,:),eye(dim))*kmrgd0';
    }];
string_jac = [
    string_sys_dimensions;
    {
    ['delays = ' delays_string ';'];
    'if any(delays < 0)';
    '    % the equation is advanced instead of delayed: we need to stop';
    '    jac = NaN(length(kmrgd));';
    '    return';
    'else';
    '    tau_max = max(delays);';
    '    if tau_max == 0';
    '        % the equation is not delayed: we need to stop';
    '        jac = NaN(length(kmrgd));';
    '        return';
    '    end';
    'end';
    };
    string_sys_nodes_etc;
    {
    '% for more details on kmrgd*, see fun_eval';
    '';
    '% compute the Jacobian matrix using finite differences for the';
    '% nonlinear part of the approximating ordinary equation (which';
    '% is defined by the right-hand side of the delay equation) and';
    '% the exact Jacobian matrix for the linear part (which is defined';
    '% by the differentiation matrix)';
    '';
    'global cds;';
    'if isfield(cds,''options'') && isfield(cds.options,''Increment'') && ~isempty(cds.options.Increment)';
    '    increment = cds.options.Increment;';
    'else';
    '    increment = 1e-5;';
    'end';
    '';
    [c_d_or_c 'jac_F_d = NaN(dim_d,length(kmrgd));'];
    [c_r_or_c 'jac_F_r = NaN(dim_r,length(kmrgd));'];
    '';
    'for ii = 1:length(kmrgd)';
    '    orig_kmrgd_ii = kmrgd(ii);';
    '    kmrgd(ii) = orig_kmrgd_ii-increment;';
    '';
    };
    string_sys_kmrgd_r_der;
    string_sys_temp;
    {
    '    % differential part of the right-hand side:';
    ['    ' c_d_or_c 'F_d1 = ['];
    };
    string_sys_eq(ind_d,1);
    {
    ['    ' c_d_or_c '];'];
    '    % renewal part of the right-hand side:';
    ['    ' c_r_or_c 'F_r1 = ['];
    };
    string_sys_eq(ind_r,1);
    {
    ['    ' c_r_or_c '];'];
    '';
    '    kmrgd(ii) = orig_kmrgd_ii+increment;';
    '';
    };
    string_sys_kmrgd_r_der;
    string_sys_temp;
    {
    '    % differential part of the right-hand side:';
    ['    ' c_d_or_c 'F_d2 = ['];
    };
    string_sys_eq(ind_d,1);
    {
    ['    ' c_d_or_c '];'];
    '    % renewal part of the right-hand side:';
    ['    ' c_r_or_c 'F_r2 = ['];
    };
    string_sys_eq(ind_r,1);
    {
    ['    ' c_r_or_c '];'];
    '';
    ['    ' c_d_or_c 'jac_F_d(:,ii) = F_d2-F_d1;'];
    ['    ' c_r_or_c 'jac_F_r(:,ii) = F_r2-F_r1;'];
    '';
    '    kmrgd(ii) = orig_kmrgd_ii;';
    'end';
    '';
    [c_d_or_c 'jac_F_d = jac_F_d/(2*increment);'];
    [c_r_or_c 'jac_F_r = jac_F_r/(2*increment);'];
    '';
    '% jac_F_r extended vertically with 0 for the differential part:';
    '% (ind_r_all-dim_d removes the offset in the indices due to the first block in kmrgd)';
    [c_r 'jac_F_r0 = repmat(jac_F_r,M,1); % only renewal equations'];
    [c_c 'jac_F_r0 = zeros(M*dim,length(kmrgd)); jac_F_r0(ind_r_all-dim_d,:) = repmat(jac_F_r,M,1); % coupled differential and renewal equations'];
    '';
    '% use sparse matrices with kron(...,eye) for performance';
    'Dkron = sparse([],[],[],M*dim,dim_d+M*dim,M*(dim_d+M*dim));';
    'speyedim = speye(dim);';
    'Dkron(:,(dim_d+1):(dim_d+M*dim)) = kron(DM(2:end,2:end), speyedim);';
    '% the first block-column is related only to the differential part';
    [c_d_or_c 'speyedim0 = speyedim(:,ind_d);'];
    [c_d_or_c 'Dkron(:,1:dim_d) = kron(DM(2:end,1), speyedim0);'];
    '';
    'jac = [';
    [c_d_or_c 'jac_F_d;'];
    [c_d 'full(Dkron); % only differential equations'];
    [c_r_or_c 'full(Dkron)-jac_F_r0; % only/also renewal equations'];
    '];';
    }];
pa_array = split(pa, ',');
string_par_values = [];
for i = 1:length(pa_array)
    string_par_values = [
        string_par_values;
        {
        sprintf('%s = parameters(%d);', pa_array{i}, i);
        }];
end
if dim_r > 0
    string_sys_rhs_r = [{
        '% --------------------------------------------------------------------------';
        'function F_r = rhs_r_eval(x,parameters,active,component,ntst,ncol)';
        '% Compute the right-hand side of renewal equations, reconstructing';
        '% the value of the corresponding components of the solution.';
        '%';
        '% x: selected columns of the x variable produced by cont; each contains';
        '%    the solution of the discretizing ODE (which may be a nontrivial orbit),';
        '%    possibly the period, and the values of the active parameters';
        '% parameters: vector of the values of the parameters on this branch';
        '%             (used for the fixed parameters)';
        '% active: indices of the active parameters';
        '% component: index (among all components) of the requested component of the';
        '%            solution of the discretizing ODE';
        '% ntst, ncol: grid parameters of the solution (ntst: number of pieces,'
        '%             ncol: number of collocation nodes in each piece)';
        '';
        'npoints = size(x,2);';
        'length_orbit = ntst*ncol+1;';
        'F_r = NaN(length_orbit,npoints);';
        '';
        'for i = 1:npoints';
        '';
        };
        string_sys_dimensions;
        {
        'dim_kmrgd = dim_d*(M+1)+dim_r*M;';
        'dim_orbit = dim_kmrgd*length_orbit;';
        'orbit = x(1:dim_orbit,i);';
        'if ntst > 0';
        '    orbit = reshape(orbit,dim_kmrgd,[]);'
        'end';
        '';
        'parameters(active) = x((end-length(active)+1):end,i);';
        };
        string_par_values;
        {
        '';
        ['delays = ' delays_string ';'];
        'if any(delays < 0)';
        '    % the equation is advanced instead of delayed: we need to stop';
        '    F_r(:) = NaN;';
        '    return';
        'else';
        '    tau_max = max(delays);';
        '    if tau_max == 0';
        '        % the equation is not delayed: we need to stop';
        '        dydt = NaN(size(kmrgd));';
        '        return';
        '    end';
        'end';
        };
        string_sys_nodes_etc;
        {
        'for j = 1:length_orbit';
        'kmrgd = orbit(:,j);';
        '';
        '% for more details on kmrgd*, see fun_eval';
        '';
        }
        string_sys_kmrgd_r_der;
        string_sys_temp;
        {
        '% renewal part of the right-hand side:';
        'switch component';
        }];
    for i = 1:dim_r
        string_sys_rhs_r = [
            string_sys_rhs_r;
            {
            sprintf('    case %d', ind_r(i));
            sprintf('        F_r(j,i) = %s', string_sys_eq{ind_r(i),1});
            }];
    end
    string_sys_rhs_r = [
        string_sys_rhs_r;
        {
        '    otherwise';
        '        F_r(j,i) = NaN;';
        'end';
        'end';
        'end';
        }];
end
string_make_initial_point_from_func = [{
    '% --------------------------------------------------------------------------';
    'function [kmrgd, tau_max] = make_initial_point_from_func(state_funcs, parameters, current_point)';
    '% Construct the initial point (vector) from the specified initial function.';
    '%';
    '% state_funcs: cell array of initial functions, in the order of the';
    '%              coordinates; each is a function of one variable in the';
    '%              delay interval';
    '% parameters: vector of the values of the parameters (possibly needed to';
    '%             compute the delays)';
    '% current_point: current value of the initial point (vector)';
    '';
    };
    string_sys_dimensions;
    {
    'kmrgd = current_point(:); % ensure it is a column vector';
    '';
    }
    string_par_values;
    {
    '';
    ['delays = ' delays_string ';'];
    'if any(delays < 0)';
    '    % the equation is advanced instead of delayed: we need to stop';
    '    kmrgd(:) = NaN;';
    '    tau_max = NaN;';
    '    return';
    'else';
    '    tau_max = max(delays);';
    '    if tau_max == 0';
    '        % the equation is not delayed: we need to stop';
    '        kmrgd(:) = NaN;';
    '        return';
    '    end';
    'end';
    };
    string_sys_nodes_etc;
    {
    'for i = 1:dim_d';
    '    if ~isempty(state_funcs{ind_d(i)})'
    '        for j = 1:M+1';
    '            kmrgd(ind_d_all(i+(j-1)*dim_d)) = state_funcs{ind_d(i)}(nodes(j));';
    '        end';
    '    end';
    'end';
    'for i = 1:dim_r';
    '    if ~isempty(state_funcs{ind_r(i)})'
    '        for j = 1:M';
    '            kmrgd(ind_r_all(i+(j-1)*dim_r)) = state_funcs{ind_r(i)}(nodes(j+1));';
    '        end';
    '        kmrgd(ind_r_all(i:dim_r:end)) = DM(2:end,2:end) \ kmrgd(ind_r_all(i:dim_r:end));';
    '    end';
    'end';
    }];
string_barint = {
    'function ff = barint(x, w, f, xx)';
    '%BARINT Barycentric interpolation';
    '%  ff = BARINT(x, w, f, xx) computes the values ff of a function at xx';
    '%  using the barycentric interpolation formula with x interpolation';
    '%  nodes, w barycentric weights and f values of the function at x.';
    '%';
    '%  Reference:';
    '%    J.-P. Berrut and L. N. Trefethen,';
    '%    Barycentric Lagrange interpolation,';
    '%    SIAM Review, 46(3):501-517, 2004,';
    '%    DOI: 10.1137/S0036144502417715';
    '';
    '% Copyright (c) 2004 Jean-Paul Berrut, Lloyd N. Trefethen';
    '% Even though the codes in the reference are not explicitly licensed,';
    '% considering that variations of them are included in Chebfun';
    '% (http://www.chebfun.org/), which is distributed under the 3-clause';
    '% BSD license, we believe that the authors'' intention was for their';
    '% codes to be freely used and that distributing them as part of MatCont,';
    '% licensed under the GPLv3 license, does not violate their rights';
    '% (the 3-clause BSD license is considered "compatible with the GPL" by';
    '% the Free Software Foundation, see';
    '% https://www.gnu.org/licenses/license-list.html).';
    '';
    'if isscalar(xx)';
    '    % most frequent cases';
    '    if xx == x(end)';
    '        ff = f(end);';
    '    elseif xx == x(1)';
    '        ff = f(1);';
    '    else';
    '        temp = w ./ (xx - x);';
    '        numer = dot(temp, f);';
    '        denom = sum(temp);';
    '        ff = numer / denom;';
    '        if isnan(ff) || isinf(ff)';
    '            j = find(xx == x);';
    '            if ~isempty(j)';
    '                ff = f(j);';
    '            end';
    '        end';
    '    end';
    'else';
    '    n = length(x);';
    '    numer = zeros(size(xx));';
    '    denom = zeros(size(xx));';
    '    exact = zeros(size(xx));';
    '    for j = 1:n';
    '        xdiff = xx - x(j);';
    '        temp = w(j) ./ xdiff;';
    '        numer = numer + temp * f(j);';
    '        denom = denom + temp;';
    '        exact(xdiff == 0) = j;';
    '    end';
    '    jj = find(exact);';
    '    ff = numer ./ denom;';
    '    ff(jj) = f(exact(jj));';
    'end';
    };
string_DE_int = {
    'function value = DE_int(f, a, b)';
    '%DE_INT Integral';
    '%  value = DE_INT(f, a, b) computes the integral of f in [a, b]';
    sprintf('%%  using the Clenshaw-Curtis quadrature formula of degree M = %d', gds.delay.QuadratureDegree);
    '%';
    '%  Reference:';
    '%    L. N. Trefethen,';
    '%    Spectral Methods in MATLAB,';
    '%    SIAM, 2000,';
    '%    DOI: 10.1137/1.9780898719598';
    '';
    '% Copyright (c) 2000 Lloyd N. Trefethen';
    '% Even though the codes in the reference are not explicitly licensed,';
    '% considering that they are made available for download on the author''s';
    '% web page and that variations of them are included in Chebfun';
    '% (http://www.chebfun.org/), which is distributed under the 3-clause';
    '% BSD license, we believe that the author''s intention was for his';
    '% codes to be freely used and that distributing them as part of MatCont,';
    '% licensed under the GPLv3 license, does not violate their rights';
    '% (the 3-clause BSD license is considered "compatible with the GPL" by';
    '% the Free Software Foundation, see';
    '% https://www.gnu.org/licenses/license-list.html).';
    '';
    '% Clenshaw-Curtis quadrature nodes and weights in [-1, 1]:';
    ['nodes = [' print_matrix_oneline_string(unit_ccnodes) '];'];
    ['weights = [' print_matrix_oneline_string(unit_ccweights) '];'];
    '% nodes rescaled to [a, b]:';
    'nodes = (b-a)/2*nodes+(a+b)/2;';
    'n = length(nodes);';
    'ff = NaN(n, 1);';
    'for i = 1:length(nodes)';
    '    ff(i) = f(nodes(i));';
    'end';
    '% scaling constant for weights: (b-a)/2';
    'value = (b-a)/2*(weights*ff);';
    };
string_aux_fun = [
    {'% --------------------------------------------------------------------------';};
    string_barint;
    {'% --------------------------------------------------------------------------';};
    string_DE_int];

fwrite=strcat(gds.system,'.m');
fwrite=fullfile(path_sys,fwrite);
[fid_write,message]=fopen(fwrite,'w');
if fid_write==-1
    errordlg(message,'Error (1)');
    return
end
fread=fullfile(path_sys,'standard.m');

[fid_read,message]=fopen(fread,'r');
if fid_read==-1
    errordlg(message,'Error (2)');
    return
end
global waithndl;
waithndl=waitbar(0,'Precomputing model');
string_handles={'out{1} = @init;';
                'out{2} = @fun_eval;';
                'out{3} = @jacobian;';
                'out{4} = [];';
                'out{5} = [];';
                'out{6} = [];';
                'out{7} = [];';
                'out{8} = [];';
                'out{9} = [];';
                'return;';};
string_init=cellstr(make_init);



waitbar(0.1);
waitbar(0.2);
waitbar(0.3);
waitbar(0.4);
waitbar(0.5);
waitbar(0.6);
waitbar(0.7);
waitbar(0.8);
waitbar(0.9);
if ~isempty(gds.options.UserfunctionsInfo)
    siz = size(gds.options.UserfunctionsInfo,2);
    for i = 1:siz
        string_handles{9+i,1}= sprintf('out{%d}= @%s;',9+i,gds.options.UserfunctionsInfo(i).name);
    end
else siz=0;end

h=0;
filecontent = '';
while feof(fid_read)==0
    tline=fgetl(fid_read);
    h=h+1;
    if h==2
        for i=1:9+siz
            fprintf(fid_write,'%s\n',string_handles{i,1});
            %filecontent = [filecontent,  sprintf('%s\n',string_handles{i,1})];
        end
        fprintf(fid_write,'delay{1} = @make_initial_point_from_func;\n');
        if dim_r > 0
            fprintf(fid_write,'delay{2} = @rhs_r_eval;\n');
        end
    end        
    matches=strrep(tline,'time,',t);
    matches=strrep(matches,'function out =','function [out, delay] =');
    matches=strrep(matches,'odefile',gds.system);
    matches=strrep(matches,',parameters',par);
    fprintf(fid_write,'%s\n',matches); 
    filecontent = [filecontent,  sprintf('%s\n',matches)]; 
    if isfield(gds,'userfunction')
        if ~isempty(findstr(matches,'varargout{1}=der5(coordinates,'))          
            for i = 1:size(gds.userfunction,2)
                hs1 = sprintf('case ''%s''\n\tvarargout{1}=%s(coordinates%s);',gds.options.UserfunctionsInfo(i).name,gds.options.UserfunctionsInfo(i).name,par);
                fprintf(fid_write,'%s\n',hs1); 
                filecontent = [filecontent,  sprintf('%s\n',hs1)]; 
            end
        end
    end        
    if ~isempty(findstr(matches,'function dydt'))
        [dim,x]=size(string_sys_dydt);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_sys_dydt{i}); 
              filecontent = [filecontent,  sprintf('%s\n',string_sys_dydt{i})];
        end
    end
    if ~isempty(findstr(matches,'handles'))
        [dim,x]=size(string_init);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_init{i});
              filecontent = [filecontent,  sprintf('%s\n',string_init{i})];
        end
    end
    if (~isempty(findstr(matches,'function jac '))&& ~isempty(string_jac))
        [dim,x]=size(string_jac);         
        for i=1:dim
              fprintf(fid_write,'%s\n',string_jac{i});
              filecontent = [filecontent,  sprintf('%s\n',string_jac{i})];
        end
    end
   
    if (~isempty(findstr(matches,'function jacp'))&& ~isempty(string_jacp))
       [dim,x]=size(string_jacp);       
       for i=1:dim
           fprintf(fid_write,'%s\n',string_jacp{i});
           filecontent = [filecontent,  sprintf('%s\n',string_jacp{i})];
       end
    end
    
    if (~isempty(findstr(matches,'function hess '))&& ~isempty(string_hess))
        [dim,x]=size(string_hess);        
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hess{i});
              filecontent = [filecontent,  sprintf('%s\n',string_hess{i})];
        end
    end
    if (~isempty(findstr(matches,'function hessp'))&& ~isempty(string_hessp))
        [dim,x]=size(string_hessp);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hessp{i});
              filecontent = [filecontent,  sprintf('%s\n',string_hessp{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens3'))&& ~isempty(string_tensor3))
        [dim,x]=size(string_tensor3);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor3{i});
              filecontent = [filecontent,  sprintf('%s\n',string_tensor3{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens4'))&& ~isempty(string_tensor4))
        [dim,x]=size(string_tensor4);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor4{i});
              filecontent = [filecontent, sprintf('%s\n',string_tensor4{i})];
        end
    end
    if (~isempty(findstr(matches,'function tens5'))&& ~isempty(string_tensor5))
        [dim,x]=size(string_tensor5);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor5{i});
              filecontent = [filecontent, sprintf('%s\n',string_tensor5{i})];
        end
    end
end
[dim,x]=size(string_make_initial_point_from_func);
for i=1:dim
    fprintf(fid_write,'%s\n',string_make_initial_point_from_func{i});
    filecontent = [filecontent,  sprintf('%s\n',string_make_initial_point_from_func{i})];
end
if dim_r > 0
    [dim,x]=size(string_sys_rhs_r);
    for i=1:dim
        fprintf(fid_write,'%s\n',string_sys_rhs_r{i});
        filecontent = [filecontent,  sprintf('%s\n',string_sys_rhs_r{i})];
    end
end
[dim,x]=size(string_aux_fun);
for i=1:dim
    fprintf(fid_write,'%s\n',string_aux_fun{i});
    filecontent = [filecontent,  sprintf('%s\n',string_aux_fun{i})];
end
newlines = strfind(filecontent, 10);
newline = newlines(1);
gds.filecontent = filecontent(newline+1:end);

if ~isempty(gds.options.UserfunctionsInfo)    
   for i=1:size(gds.options.UserfunctionsInfo,2)
       str_user = []; res=0;
       if isfield(gds,'userfunction') && ~isempty(gds.userfunction{i})
           str_user = import_delay('replace_delays', cellstr(renameforsym(gds.userfunction{i}, original_cor, original_pa)), num_delays);
       else 
           str_user=cellstr('res=');
       end
       hs1 = sprintf('function userfun%d=%s(t,kmrgd%s)',i,gds.options.UserfunctionsInfo(i).name,par);
       fprintf(fid_write,'%s\n',hs1);
       hs1 = sprintf('userfun%d',i);
       dim = size(str_user,1);
       for j = 1:dim
           userline = str_user{j};
           d = strmatch('res=',userline,'exact');
           if findstr('res',userline),res=1;end
           userline = strrep(userline,'res',hs1);
           if d==1
               fprintf(fid_write,'\t%s=0;\n',hs1);
           else 
               fprintf(fid_write,'\t%s;\n',userline);
           end
       end
       if res==0,fprintf(fid_write,'\t%s=0;\n',hs1);end
   end
end            
    
waitbar(0.95);
fclose(fid_read);
fclose(fid_write);
file=fullfile(path_sys,gds.system);

%fix mapping here
gds.delay.equations = original_equations;
gds.delay.coordinates = [toGdsStruct(original_cor), gds.delay.coordinates(:,3)];
gds.equations = char({
    '% This system was created with the delay equations importer.';
    '% DO NOT EDIT THIS SYSTEM with MatCont''s standard system editor.';
    '% Use the delay equations editor instead.';
    '% The following are the equations originally input in the delay equations importer.';
    gds.delay.equations});
gds.coordinates = create_aux_variables;
gds.parameters = toGdsStruct(original_pa);
gds.class = 'delay equation (imported as ODE)';

save(file,'gds');
gds.ok = true;
hds=[];
delete(waithndl);
delete(handles.system);
%set(MC.mainwindow.compute,'enable','off');
%set(MC.mainwindow.window,'enable','off');
%set(MC.mainwindow.Type,'enable','on');
%set(MC.mainwindow.select_userfunctions,'enable','on');
%{
if gds.open.figuur==1;starter;end
if gds.open.continuer==1;continuer;end
if gds.open.numeric_fig==1;numeric;end
if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
if gds.open.integrator==1;integrator;end
%}
cd(path_sys);cd ..;
rehash;

% remove old multilinear if exist
multilinearforms_file = ['./Systems/', gds.system, '_multilinearforms.m'];
if exist(multilinearforms_file, 'file')
    delete(multilinearforms_file); 
end


% tempstr.label = 'Point';
% tempstr.Tag = 'P_O_DO';
% matcont('point_callback',tempstr)

% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel.
global gds oldgds path_sys;
gds=oldgds;
if (~isempty(gds))
    file=fullfile(path_sys,gds.system);
    save(file,'gds');
    load_system(handles);
    %if gds.open.figuur==1;starter;end
    %if gds.open.continuer==1;continuer;end
    %if gds.open.numeric_fig==1;numeric;end
    %if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
    %if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
    %if gds.open.integrator==1;integrator;end
end
delete(handles.system);

% --------------------------------------------------------------------
function varargout = coordinates_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.coordinates.
global gds;
% store old coordinates, they might contain usefull initial values
old_coord = gds.delay.coordinates;
gds.delay.coordinates=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC='';
str=get(handles.coordinates,'String');
str=parse(str);
t=1;
for i=1:length(str)
    co{t,1}=str{i,1};
    string=str{i};
    str1=findstr(string,'[');
    str2=findstr(string,']');
    m=length(str1);
    if (m==length(str2)&&(m==1))
        num=str2double(string(str1+1:str2-1));
        var=string(1:str1-1);
        co{t,1}=strcat(var,'(1)');
        if num>1
            for j=2:num
                t=t+1;
                co{t,1}=strcat(var,'(',num2str(j),')');
            end
        end
    end
    t=t+1;
end
gds.delay.coordinates=co;
gds.delay.dim=size(gds.delay.coordinates,1);
if ((gds.delay.dim==1)&&(strcmp(gds.delay.coordinates{1},'')))
    gds.delay.coordinates=[];
    gds.delay.dim=0;
else
    % run through each coordinate and try to set its initial value
    % from the old coordinate
    for i=1:gds.delay.dim
        gds.delay.coordinates{i,2}=0;
        for j = 1:size(old_coord,1)
            if strcmp(old_coord{j,1},gds.delay.coordinates{i,1})
                gds.delay.coordinates{i,2}=old_coord{j,2};
                break;
            end
        end
    end
end




% --------------------------------------------------------------------
function varargout = parameters_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.parameters.
global gds;
% store old parameters, they might contain usefull initial values
old_params = gds.parameters;
gds.parameters=[];gds.options.ActiveParams=[];
gds.T=[];gds.eps0=[];gds.eps1=[];gds.extravec=[1 0 0]; gds.t=[]; gds.epsilon=[];
gds.options.IgnoreSingularity=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC='';
str=get(handles.parameters,'String');
gds.parameters=parse(str);
[ndim,x]=size(gds.parameters);


if ((ndim==1)&&(strcmp(gds.parameters{1},'')))
    gds.parameters=[];
else
    % run through each parameter and try to set its initial value
    % from the old parameter
    for i=1:ndim
        gds.parameters{i,2}=0;
        for j = 1:size(old_params,1)
            if strcmp(old_params{j,1},gds.parameters{i,1})
                gds.parameters{i,2}=old_params{j,2};
                break;
            end
        end
    end
end

  
% --------------------------------------------------------------------
function varargout = time_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.time.
global gds
gds.time=[];
str=get(handles.time,'String');
gds.time={str, 0};

% --------------------------------------------------------------------
function varargout = name_system_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.name_system.
global gds path_sys
str1=get(handles.name_system,'String');
str1=strrep(str1,'.m','');
str=strcat(str1,'.m');
file=fullfile(path_sys,str);
if (exist(file)==2)
    warning_box(str1,file,handles);
else 
    w=which(str);
    if ~isempty(w)
        errordlg('This name is already in use for another function (Possibly in Matlab itself).','wrong name');
        set(handles.name_system,'String','');gds.system='';
        return
    end
    gds.system=str1;
end

%--------------------------------------------------------------------
function varargout=load_system(handles)
global gds
co='';par='';
if ((gds.delay.dim)~=0)
    co=gds.delay.coordinates{1,1};
end
if ((gds.delay.dim)>=2)
  for i=2:gds.delay.dim
       co=horzcat(co,',',gds.delay.coordinates{i,1});
  end
end
dim=size(gds.parameters,1);
if (dim~=0)
    par=gds.parameters{1,1};
end
if (dim>=2)
   for i=2:dim
        par=horzcat(par,',',gds.parameters{i,1});
   end
end  
set(handles.name_system,'String',gds.system);
set(handles.coordinates,'String',co);
set(handles.parameters,'String',par);
set(handles.time,'String',gds.time{1,1});
% before setting the system, trim the equations
equations=gds.delay.equations;
equations_string = '';
for equation = equations'
    equations_string = sprintf('%s%s\n',equations_string, strtrim(equation'));
end
% remove last newline character
equations_string = equations_string(1:end-1);
set(handles.sys,'String',equations_string);
set(handles.collocation_degree,'String',gds.delay.CollocationDegree);
set(handles.quadrature_degree,'String',gds.delay.QuadratureDegree);
guidata(handles.system,handles);

%---------------------------------------------------------------------
function all=parse(input_string)
global gds;
remainder=input_string;
all='';
while (any(remainder))
[chopped,remainder]=strtok(remainder,',');
all=strvcat(all,chopped);
end
all=cellstr(all);

%--------------------------------------------------------------------
function [string_temp, string_eq]=replace_sys_input(num_delays)
global gds;
string = cellstr(gds.delay.equations);
if isempty(string)
    return
end
string_temp = '';
string_eq = '';
[temp,eq] = parse_input(string);
num_temp = size(temp,1);
for i = 1:num_temp
    string_temp{i,1} = strcat(temp{i,1},';');
end
num_eq = size(eq,1);
for i=1:num_eq
    string_eq{i,1} = '';
end
for j = 1:num_eq
      [~,string_eq{j,1}] = strtok(eq{j,1},'=');
      string_eq{j,1} = strtok(string_eq{j,1},'=');
      string_eq{j,1} = strcat(string_eq{j,1},';');
end
if (num_eq~=gds.delay.dim)
  error('The left-hand sides do not match the coordinates.');
end
string_temp = replace_delays(string_temp, num_delays);
string_eq = replace_delays(string_eq, num_delays);
%---------------------------------------------------------------------------
function [temp,eq] = parse_input(string)
global gds;
% can't parse empty strings, so first remove them here
j=1;
cleaned_string={};
for i=1:size(string,1)
    if ~strcmp(string{i},'')
        cleaned_string{j}=string{i};
        j=j+1;
    end
end
string=cellstr(cleaned_string');
% continue with code
dim = size(string,1);
vars = size(gds.delay.coordinates,1);
temp='';eq='';
p=1;s=1;
for j=1:dim
    if (~contains(string{j},'='))
        error('All lines must be either assignments to temporary variables or equations.');
    end
    k=[];  
    for i=1:vars
        teststring = string{j};
        if exist('strtrim','builtin')
            coordinate = strtrim(gds.delay.coordinates{i,1});
        else
            coordinate = deblank(gds.delay.coordinates{i,1});
        end
        match = strcat('\<',coordinate,'\>');
        [pos, endpos] = regexp(teststring,match);
        if ~isempty(pos) && pos(1)==1
           k = 1;
           if (vars-i ~= dim-j)
               error('Equations are in the wrong order, compared to the coordinates.');
           end
           nextchar = teststring(endpos(1)+1);
        end
    end
    if isempty(k)
        temp{p,1} = string{j};
        p = p+1;
    else
        eq{s,1} = string{j};
        if (nextchar=='''')
            % differential equation (ordinary or delayed)
            gds.delay.coordinates{s,3} = 1;
        else
            % renewal equation
            gds.delay.coordinates{s,3} = 0;
        end
        s = s+1;
    end
end

%------------------------------------------------------------------------
function [num_delays, delays_string] = extract_delays
global gds;
delays = [];
for ind_eq=1:size(gds.delay.equations,1)
    % get right-hand side of current equation
    [rhs,~] = split(gds.delay.equations(ind_eq,:),'=');
    rhs = rhs{2};
    % syntax for distributed delays:
    % $\int_a^b g(x(t+\theta)) d\theta$ -> DE_int(@(theta)g(x[t+theta]), a, b)
    expression = 'DE_int\(((?!\<DE_int\>).)*,((?!\<DE_int\>).)*,((?!\<DE_int\>).)*\)';
    matches = regexp(rhs,expression,'match');
    matches = unique(matches);
    for ind_mat=1:length(matches)
        % correctly match opening and closing parentheses of DE_int
        ind_open_close = find_matching_parentheses(matches{ind_mat});
        [ind_open_row, ~] = find(ind_open_close==7);
        arguments_string = extractBetween(matches{ind_mat},ind_open_close(ind_open_row,1),ind_open_close(ind_open_row,2));
        % find endpoints of distributed delay interval
        integration_variable = regexp(arguments_string{1},'@\(([^\)]*)\)','tokens');
        delays_in_integrand = regexp(arguments_string{1},['\[' gds.time{1,1} '([^\]\[]*' integration_variable{1}{1} '[^\]\[]*)\]'],'tokens');
        tokens = regexp(arguments_string{1},'\((.*),([^,]*),([^,]*)\)','tokens');
        new_arguments_string = arguments_string{1};
        for ind_del=1:length(delays_in_integrand)
            delays = [delays {strrep(delays_in_integrand{ind_del}{1},integration_variable{1}{1},['(' strtrim(tokens{1}{2}) ')'])}];
            delays = [delays {strrep(delays_in_integrand{ind_del}{1},integration_variable{1}{1},['(' strtrim(tokens{1}{3}) ')'])}];
            new_arguments_string = strrep(new_arguments_string,['[' gds.time{1,1} delays_in_integrand{ind_del}{1} ']'],'');
        end
        % delete argument of distributed delays to avoid spurious discrete delays
        rhs = strrep(rhs,['DE_int' arguments_string{1}],['DE_int' new_arguments_string]);
    end
    % syntax for discrete delays: x[t-1], x[t-a], x[t+b]
    for ind_cor=1:size(gds.delay.coordinates,1)
        % x[t-...] -> /x\[t(\W[^\]\[]*)\]/ -> '-...'
        expression = [gds.delay.coordinates{ind_cor,1} '\[' gds.time{1,1} '(\W[^\]\[]*)\]'];
        tokens = regexp(rhs,expression,'tokens');
        delays = [delays tokens{:}];
    end
end
delays = unique(delays);
num_delays = length(delays);
if num_delays == 0
    error('No delay found.')
end
start_minus = startsWith(delays,'-');
start_plus = startsWith(delays,'+');
for ind_del = 1:num_delays
    if start_minus(ind_del)
        delays{ind_del}(1) = '';
    elseif start_plus(ind_del)
        delays{ind_del}(1) = '-';
    else
        error('Delays must be specified in the form [t-expr] or [t+expr] (where t is the time variable).')
    end
end
delays_string = ['[' delays{1}];
for ind_del=2:num_delays
    delays_string = [delays_string ',' delays{ind_del}];
end
delays_string = [delays_string ']'];


%------------------------------------------------------------------------
function check_delays(delays_string)
global gds;
for i=1:size(gds.parameters,1)
    eval(sprintf('%s=%d;\n',gds.parameters{i,1},gds.parameters{i,2}));
end
clear i gds
eval([delays_string, ';']);
%------------------------------------------------------------------------
function pairs = find_matching_parentheses(string)
% pairs = find_matching_parentheses(string)
% finds pairs of matching parentheses in string and returns the
% corresponding pairs of indices. Parentheses are matched starting from the
% innermost and unmatched parentheses are discarded.

% Function by Daniel Renjewski, 15 March 2023, retrieved from
% https://it.mathworks.com/matlabcentral/answers/121920-how-do-i-match-nested-parenthesis-brackets-or-braces-with-dynamic-regular-expressions
% (license: Creative Commons Attribution-ShareAlike 3.0 Unported,
% CC-BY-SA 3.0) adapted to tolerate and discard unmatched parentheses
% (license of the adaptation: Creative Commons Attribution-ShareAlike 4.0
% International, CC-BY-SA 4.0; the license upgrade is allowed by section
% 4(b) of CC-BY-SA 3.0). The inclusion of the adapted version in this
% program, licensed under the GNU General Public License version 3 (GPLv3),
% is allowed since the GPLv3 has been declared by the Creative Commons as a
% "BY-SA–Compatible License" (see section 3(b) of CC-BY-SA 4.0 and
% https://creativecommons.org/compatiblelicenses).

pairs = [];
opening = strfind(string,'(');
closing = strfind(string,')');
while ~(isempty(opening)||isempty(closing))
    matching = find(opening<closing(1),1,'last');
    if ~isempty(matching)
        pairs = [pairs; opening(matching), closing(1)];
        opening(matching) = [];
    end
    closing(1) = [];
end

%------------------------------------------------------------------------
function strings = replace_delays(strings,num_delays)
global gds;
[ind_d,ind_r] = create_indices;
% same treatment for discrete delays and coordinates in distributed delays
for ind_str=1:length(strings)
    string = strings{ind_str};
    % syntax for delays: x[t-1], x[t-a], x[t+b], x[t+theta]
    for ind_cor=1:size(gds.delay.coordinates,1)
        cor = gds.delay.coordinates{ind_cor,1};
        cor_type = gds.delay.coordinates{ind_cor,3};
        if cor_type == 1 % differential equation (ordinary or delayed)
            ind_cor0 = find(ind_d==ind_cor,1); % index in the current-time part of kmrgd
            % replace current time terms with explicit time dependency...
            string = strrep(string,[cor '[' gds.time{1,1} ']'],sprintf('kmrgd(%d)',ind_cor0));
            % ... and without it
            string = regexprep(string,[cor '($|[^\[a-zA-Z0-9_])'],sprintf('kmrgd(%d)$1',ind_cor0));
            % x[t-...] -> /x\[t(\W[^\]\[]*)\]/ -> '-...'
            expression = [cor '\[' gds.time{1,1} '(\W[^\]\[]*)\]'];
            matches = regexp(string,expression,'match');
            matches = unique(matches);
            for ind_mat=1:length(matches)
                if num_delays == 1
                    % every delay must be tau_max: we can use the last block of kmrgd
                    string = strrep(string,matches{ind_mat},sprintf('kmrgd(dim_d+(M-1)*dim+%d)',ind_cor));
                else
                    tokens = regexp(matches{ind_mat},expression,'tokens');
                    string = strrep(string,matches{ind_mat},sprintf('barint(nodes,weights,kmrgd([%d,(dim_d+%d):dim:(dim_d+M*dim)]),%s)',ind_cor0,ind_cor,tokens{1}{1}));
                end
            end
        elseif cor_type == 0 % renewal equation
            % replace also discrete delays, although that does not make
            % much sense and poses numerical problems (derivative...)
            ind_cor0 = find(ind_r==ind_cor,1); % index in the current-time part of kmrgd_r_der
            % replace current time terms with explicit time dependency...
            string = strrep(string,[cor '[' gds.time{1,1} ']'],sprintf('kmrgd_r_der(%d)',ind_cor0));
            % ... and without it
            string = regexprep(string,[cor '($|[^\[a-zA-Z0-9_])'],sprintf('kmrgd_r_der(%d)$1',ind_cor0));
            % x[t-...] -> /x\[t(\W[^\]\[]*)\]/ -> '-...'
            expression = [cor '\[' gds.time{1,1} '(\W[^\]\[]*)\]'];
            matches = regexp(string,expression,'match');
            matches = unique(matches);
            for ind_mat=1:length(matches)
                if num_delays == 1
                    % every delay must be tau_max: we can use the last block of kmrgd_r_der
                    string = strrep(string,matches{ind_mat},sprintf('kmrgd_r_der(M*dim_r+%d)',ind_cor0));
                else
                    tokens = regexp(matches{ind_mat},expression,'tokens');
                    string = strrep(string,matches{ind_mat},sprintf('barint(nodes,weights,kmrgd_r_der(%d:dim_r:((M+1)*dim_r)),%s)',ind_cor0,tokens{1}{1}));
                end
            end
        end
    end
    strings{ind_str} = string;
end

%------------------------------------------------------------------------
function [ind_d,ind_r,dim_d,dim_r,ind_d_all,ind_r_all] = create_indices
% ind_{d,r}: indices of coordinates described by {differential,renewal} equations
% dim_{d,r}: dimension of the {differential,renewal} part of the system
% ind_{d,r}_all: indices of all components of kmrgd described by {differential,renewal} equations
global gds;
ind_d = find(cell2mat(gds.delay.coordinates(:,3))==1);
ind_r = find(cell2mat(gds.delay.coordinates(:,3))==0);
dim_d = length(ind_d);
dim_r = length(ind_r);
ind_d_all = [1:dim_d, NaN(1,dim_d*gds.delay.CollocationDegree)];
ind_r_all = NaN(1,dim_r*gds.delay.CollocationDegree);
for i=1:gds.delay.CollocationDegree
    ind_d_all((1:dim_d)+i*dim_d) = dim_d+ind_d+(i-1)*gds.delay.dim;
    ind_r_all((1:dim_r)+(i-1)*dim_r) = dim_d+ind_r+(i-1)*gds.delay.dim;
end

%------------------------------------------------------------------------
function vars = create_aux_variables
% order of variables in the state:
% [differential part at current time (length dim_d);
%  full state at first node in the past (length dim);
%  ...
%  full state at M-th node in the past (length dim)]
global gds;
dim = gds.delay.dim;
[ind_d,~,dim_d] = create_indices;
M = gds.delay.CollocationDegree;
n_digits = floor(log10(M))+1;
vars = cell(dim_d+M*dim,2);
for d = 1:dim
    if gds.delay.coordinates{d,3} == 1
        d0 = find(ind_d==d,1);
        vars(d0,:) = gds.delay.coordinates(d,1:2);
    end
    for m = 1:M
        vars(dim_d+(m-1)*dim+d,:) = gds.delay.coordinates(d,1:2);
        vars{dim_d+(m-1)*dim+d,1} = sprintf(sprintf('%%s_aux%%0%dd', n_digits), vars{dim_d+(m-1)*dim+d,1}, m);
    end
end

%------------------------------------------------------------------------
function ms = print_matrix_oneline_string(m)
rows = size(m,1);
cols = size(m,2);
ms = '';
for i=1:rows
    if cols > 0
        ms = [ms sprintf('%.15g',m(i,1))];
    end
    if cols > 1
        ms = [ms sprintf(',%.15g',m(i,2:end))];
    end
    if i < rows
        ms = [ms ';'];
    end
end

%------------------------------------------------------------------------
function ms = print_matrix_multiline_cell(m)
rows = size(m,1);
cols = size(m,2);
ms = cell(rows,1);
for i=1:rows
    mss = print_matrix_oneline_string(m(i,:));
    ms{i} = [mss ';'];
end

%------------------------------------------------------------------------
function [x, w, D] = cheb(N)
%CHEB Chebyshev type II nodes (extrema) and relevant quantities
%  [x, w, D] = CHEB(N) returns the N+1 Chebyshev nodes x in [-1, 1]
%  and the corresponding barycentric interpolation weights w and
%  differentiation matrix D.
%
%  References:
%    L. N. Trefethen,
%    Spectral Methods in MATLAB,
%    SIAM, 2000,
%    DOI: 10.1137/1.9780898719598
%
%    J.-P. Berrut and L. N. Trefethen,
%    Barycentric Lagrange interpolation,
%    SIAM Review, 46(3):501-517, 2004,
%    DOI: 10.1137/S0036144502417715

assert(N>0)
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D-diag(sum(D,2));

w = [1/2; ones(N-1,1); 1/2].*(-1).^(0:N)';

%------------------------------------------------------------------------
function [w, x] = clencurt(N)
%CLENCURT Weights and nodes of Clenshaw-Curtis quadrature
%  [w, x] = CLENCURT(N) returns the weights of the Clenshaw-Curtis
%  quadrature, i.e. the pseudospectral method on N+1 Chebyshev
%  nodes x in [-1, 1].
%
%  Reference:
%    L. N. Trefethen,
%    Spectral Methods in MATLAB,
%    SIAM, 2000,
%    DOI: 10.1137/1.9780898719598

theta = pi * (0:N)' / N;
x = cos(theta);
w = zeros(1, N+1);
ii = 2:N;
v = ones(N-1, 1);
if mod(N, 2) == 0
    w(1) = 1 / (N^2-1);
    w(N+1) = w(1);
    for k = 1:N/2-1
        v = v - 2 * cos(2*k*theta(ii)) / (4*k^2-1);
    end
    v = v - cos(N*theta(ii)) / (N^2-1);
else
    w(1) = 1 / N^2;
    w(N+1) = w(1);
    for	k = 1:(N-1)/2
        v = v - 2 * cos(2*k*theta(ii)) / (4*k^2-1);
    end
end
w(ii) = 2 * v / N;

%------------------------------------------------------------------------
function init
global gds;
    gds = []; gds.delay.coordinates = []; gds.coordinates = []; gds.parameters = [];
    gds.time{1,1} = 't';gds.time{1,2} = 0; gds.options = contset;
    gds.system = '';
    gds.class = '';
    gds.curve.new = '';gds.curve.old = '';
    gds.delay.equations = []; gds.equations = [];
    gds.delay.dim = 0; gds.dim = 0;
    gds.der = [[1 1 1 1 1];zeros(3,5)]; 
    gds.jac = '';%string that contains the jacobian
    gds.jacp = '';%string that contains the jacobianp
    gds.hess = '';%string that contains the hessian
    gds.hessp = '';%string that contains the hessianp
    gds.tensor3 = ''; gds.tensor4 = ''; gds.tensor5 = '';
    gds.point = ''; gds.type = '';
    gds.discretization.ntst = 20; gds.discretization.ncol = 4;
    gds.period = 1;
    gds.plot2 = '';gds.plot3 = '';gds.PRC='';gds.dPRC='';
    gds.open.figuur = 0; gds.open.continuer = 0; gds.open.numeric_fig = 0;
    gds.open.D2 = 0;gds.open.D3 = 0;gds.open.PRC = 0; gds.open.dPRC = 0; gds.open.integrator = 0;
    gds.integrator = []; gds.integrator.method = 'ode45'; gds.integrator.options = [];
    gds.integrator.tspan = [0 1]; gds.numeric = [];
    gds.numeric.O = {'time' 1;'coordinates' 1;'parameters' 0'};
    gds.numeric.EP = {'coordinates' 1;'parameters' 1;'testfunctions' 0;'eigenvalues' 0;'current stepsize' 0};
    % XXXX
    %gds.numeric.LC = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.LC = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0;'PRC' 0;'dPRC' 0;'Input' 0};
    % XXXX
    gds.numeric.PD = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    % XXX
    gds.numeric.Hom = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.numeric.HSN = {'parameters' 1;'period' 1;'testfunctions' 0;'multipliers' 0;'current stepsize' 0};
    gds.diagram = 'diagram';
    gds.parameters=[];gds.options.ActiveParams=[];
    gds.T=[];gds.eps0=[];gds.eps1=[];gds.extravec=[1 0 0]; gds.t=[]; gds.epsilon=[];
    gds.options.IgnoreSingularity=[];gds.plot2='';gds.plot3='';gds.PRC='';gds.dPRC=''; 
    gds.options.PRC = 0; gds.options.dPRC = 0;
    gds.delay.CollocationDegree=10; gds.delay.QuadratureDegree=10;
         
%---------------------------------------------------------------------
function warning_box(stri1,stri2,handles)
global gds path_sys
button = questdlg('System already exist! Do you want to continue? If you press yes to continue, you will overwrite the existing system',...
'System already exist','Yes','No','No');
dir=path_sys;
if strcmp(button,'Yes')
   gds.system=stri1;
   delete(stri2);
elseif strcmp(button,'No')
   set(handles.name_system,'String','');
end

%---------------------------------------------------------------------------
function string=make_init
global gds;
string{1,1}='y0=[';
if (gds.dim>1)
    for i=1:(gds.dim-1)
        string{1,1}=strcat(string{1,1},'0,');
    end
end
string{1,1}=strcat(string{1,1},'0];');
string{2,1}=strcat('options = odeset(''Jacobian'',[],''JacobianP'',[],''Hessians'',[],''HessiansP'',[]);');   

%-----------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function name_system_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name_system (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function coordinates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sys_Callback(hObject, eventdata, handles)
% hObject    handle to sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sys as text
%        str2double(get(hObject,'String')) returns contents of sys as a double


% --- Executes during object creation, after setting all properties.
function sys_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- toGdsStruct('x,y,z')  ->      {'x'    [0], 'y'    [0], 'z'    [0]}
function result = toGdsStruct(str)
    items = strsplit(str, ',');
    result = [items', num2cell(zeros(length(items), 1))];




function collocation_degree_Callback(hObject, eventdata, handles)
% hObject    handle to collocation_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of collocation_degree as text
%        str2double(get(hObject,'String')) returns contents of collocation_degree as a double
global gds
val = abs(round(str2double(get(handles.collocation_degree,'String'))));
if val>0
    set(handles.collocation_degree,'String',val);
    gds.delay.CollocationDegree = val;
else
    set(handles.collocation_degree,'String',gds.delay.CollocationDegree);
end


% --- Executes during object creation, after setting all properties.
function collocation_degree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to collocation_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function quadrature_degree_Callback(hObject, eventdata, handles)
% hObject    handle to quadrature_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quadrature_degree as text
%        str2double(get(hObject,'String')) returns contents of quadrature_degree as a double
global gds
val = abs(round(str2double(get(handles.quadrature_degree,'String'))));
if val>0
    set(handles.quadrature_degree,'String',val);
    gds.delay.QuadratureDegree = val;
else
    set(handles.quadrature_degree,'String',gds.delay.QuadratureDegree);
end


% --- Executes during object creation, after setting all properties.
function quadrature_degree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quadrature_degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in info_syntax.
function info_syntax_Callback(hObject, eventdata, handles)
% hObject    handle to info_syntax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
info_syntax

