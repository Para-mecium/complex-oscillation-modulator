classdef FMAM_ODE < handle
    properties (Constant)
       maxiternum = 100; 
       maxerr = 1;
       errBound = 1e-8;
       Lconst = 500;
       coe_momentum = 0;
       coeConst = 1e3;
       compM = 20;
    end

    properties 
        sys % the functions of the dynamical system, expected to be a 1*N cell
        obs % the functions of observables, expected to be a 1*n cell
        derivatives % the derivatives of system function and observables
        items_perturb % quantities which need to be modulated, expected to be a 1*n_per struct 
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %                 Attributes                   %
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %   -- prop:      propname(str)    
                      %
                      %      proplist = {params, p_Psi, q_Psi, p_var, q_var, p_obs, 
                      %                  q_obs, varPhiMax, obsPhiMax,
                      %                  varPhiMin, obsPhiMin, varAmp,
                      %                  obsAmp, varMax, obsMax, varMin, 
                      %                  obsMin, varPhase, obsPhase}
                      %
                      %   -- idx:       index(int)
                      %
                      %      format:    p_var(q_var), p_obs(q_obs)
                      %                      -- [i,j] for jth the Fourier
                      %                      coefficient of ith term
                      %                 varPhase(obsPhase) 
                      %                      -- [i,j] for the phase difference 
                      %                         between var_i(obs_i) and 
                      %                         var_j(obs_j)
                      %                 others -- i 
                      %
                      %   -- target:    modulation target(double)
                      %   For props that have only one single component,
                      %   the 'idx' attr should be set as 1.
                      %
                      %   E.g. iterms_perturb(1,i) = struct('prop',
                      %   'varAmp', 'idx', 4, 'target', 10) means that the
                      %   ith modulation target is to vary the amplitude of
                      %   the 4th system variable to 10.

        items_controlled % quantities which is varied to achieve the modulation target, 
                         % expected to be a vector
        maxstepsize % ith component is the maxstepsize of ith quantity of items_perturb
        truncationOrder % truncation order of trigonometric polynomials
        accuracy % 
        stat % instance of the class 'state'
        maxLr = 1 % maximum learning rate
        minLr = 1e-4 % minimum learning rate
        lr % learning rate  
        isPsiUpdated = false 

        needLog = false
        logs % used to save the info of each steps
    end

    properties (Access = private)
        dimVar
        dimObs
        dimParams 
        n_per
        stepsize = 0
        x0Last
        paramsLast
        isbreak = false
        
        p_Psi_init = []
        q_Psi_init = []
    end

    properties(Dependent)
        items_per_curr
    end

    methods
        function obj = FMAM_ODE(system,observables,stat,items_per,Coe_Controlled,maxstepsize,err)
            obj.sys = system;
            obj.obs = observables;
            obj.lr = obj.maxLr;
            obj.items_perturb = items_per;
            obj.items_controlled = Coe_Controlled;
            obj.maxstepsize = maxstepsize;
            obj.accuracy = err;

            obj.stat = stat;
            obj.truncationOrder = stat.truncationOrder;
            
            obj.dimVar = size(system,2);
            obj.dimObs = size(observables,2);
            obj.dimParams = size(stat.params,2);
            obj.n_per = size(items_per,2);

            obj.x0Last = sum(obj.stat.p_var,1);
            obj.paramsLast = obj.stat.params;

            if ~obj.isPsiUpdated
                obj.p_Psi_init = obj.stat.p_Psi;
                obj.q_Psi_init = obj.stat.q_Psi;
            end

            % Compute the derivatives of the right-hand side of the ODE
            obj.derivatives = struct;
            Derivative_system = struct;
            parameters_sym = sym('parameter',[1,obj.dimParams]);
            variable_sym = sym('variable',[1,obj.dimVar]);
            
            for i = 1:obj.dimVar
                fi = system{i};
                J = jacobian(fi(variable_sym,parameters_sym),[variable_sym,parameters_sym]);
                for j = 1:obj.dimVar + obj.dimParams
                    Derivative_system(i,j).function = matlabFunction(J(1,j),'vars',{variable_sym,parameters_sym});
                end
            end
            obj.derivatives.var = Derivative_system;
            clearvars parameters_sym variable_sym J

            % Compute the derivatives of the observables
            Derivative_observable = [];
            DDerivative_observable = [];
            if obj.dimObs ~= 0
                Derivative_observable = struct;
                DDerivative_observable = struct;
                variable_sym = sym('variable',[1,obj.dimVar]);
                for i = 1:obj.dimObs
                    fi = observables{1,i};
                    J = jacobian(fi(variable_sym),variable_sym);
                    for j = 1:obj.dimVar
                        dfidxj = matlabFunction(J(1,j),'vars',{variable_sym});
                        Derivative_observable(i,j).function = dfidxj;
                        JJ = jacobian(dfidxj(variable_sym),variable_sym);
                        for k = 1:obj.dimVar
                            % d^2 fi/(dxjdxk)
                            DDerivative_observable(i,j,k).function = (matlabFunction(JJ(1,k),'vars',{variable_sym}));
                        end
                    end
                end
            end
            obj.derivatives.obs = Derivative_observable;
            obj.derivatives.obs2 = DDerivative_observable;
            clearvars variable_sym J JJ
        end

        function val = get.items_per_curr(obj)
            val  = zeros(1,obj.n_per);
            for i = 1:length(val)
                item_per = obj.items_perturb(i);
                propname = item_per.prop;

                index = item_per.idx;
                if length(index) > 1
                    index = num2cell(index);
                    index = sub2ind(size(obj.stat.(propname)),index{:});
                end
                val(i) = obj.stat.(propname)(index);
            end
        end

        function perturb(obj)
            for i = 1:obj.n_per
                item_per = obj.items_perturb(i);
                propname = item_per.prop;
                index = item_per.idx;
                if length(index) > 1
                    index = num2cell(index);
                    index = sub2ind(size(obj.stat.(propname)),index{:});
                end
                obj.stat.add(propname,index,obj.stepsize(i))
            end
        end

        function step(obj)
            obj.autostepsize();
            isdone = all(obj.stepsize == 0);
            count = 0;
            obj.isbreak = false;
            while ~isdone
                obj.lr = obj.maxLr;
                count = count+1;
                disp(['Computing ' num2str(count) 'th iteration...'] )
                % update from the latest version
                obj.x0Last = sum(obj.stat.p_var,1);
                obj.paramsLast = obj.stat.params;
                % perturbation
                obj.perturb()
                % Increment
                obj.fit();
                % obj.stat.updatePeriod();
                % obj.stat.updateVar2();
                if obj.isbreak
                    break
                end

                if obj.needLog
                    log_curr = struct('params',obj.stat.params);
    
                    for i = 1:obj.n_per
                        item_per = obj.items_perturb(i);
                        propname = item_per.prop;
        
                        index = item_per.idx;
                        if length(index) > 1
                            index = num2cell(index);
                            index = sub2ind(size(obj.stat.(propname)),index{:});
                        end
                        
                        logName = [propname '_idx_' num2str(index)];
                        log_curr.(logName) = obj.stat.(propname)(index);
                    end
    
                    obj.logs = [obj.logs log_curr];
                end

                obj.autostepsize();
                isdone = all(obj.stepsize == 0);
            end
        end

        function autostepsize(obj)
            init_value = obj.items_per_curr;
            maxstep = obj.maxstepsize;
            ExpectedError = obj.accuracy;

            assert(length(init_value) == length(maxstep));

            step = zeros(1,length(init_value));

            for i = 1:length(init_value)
                end_value_i = obj.items_perturb(i).target;
                direction = 2*(end_value_i > init_value(i))-1;
                if abs(end_value_i - init_value(i))>ExpectedError
                    exponent=floor(log(abs(end_value_i - init_value(i)))/log(10));
                    step(i) = direction*min(10^(exponent),maxstep(i));
                else
                    step(i) = 0;
                end
            end
            obj.stepsize = step;
        end

        function fit(obj)
            err = obj.res();
            count = 0;
            while any(err > FMAM_ODE.errBound) 
                count = count + 1;

                obj.oneIter();
                err = obj.res();
                disp(['Error: ' num2str(err)])
                if count > obj.maxiternum
                    disp(['Iteration number exceeds maxiternum(' num2str(obj.maxiternum) ')!'])
                    disp(['Error: ' num2str(err)])
                    obj.isbreak = true;
                    break
                elseif all([count > 1, any(err > obj.maxerr)])
                    disp(['Error exceeds maxerr(' num2str(obj.maxerr) ')!'])
                    disp(['Iteration number: ' num2str(count)])
                    obj.isbreak = true;
                    break
                end
            end
        end

        function val = res(obj)
            system = obj.sys;

            p_var = obj.stat.p_var;
            q_var = obj.stat.q_var;
            p_Psi = obj.stat.p_Psi;
            q_Psi = obj.stat.q_Psi;
            parameters = obj.stat.params;
            % x0_curr = sum(p_var,1);
            % x0_last = obj.x0Last;
            % params_last = obj.paramsLast;

            L = obj.Lconst;
            N = obj.dimVar;
            % n = obj.dimObs;
            M = obj.truncationOrder;

            phi = (0:L-1)'*2*pi/L;
            [vc,vs] = obj.Vec_CS(phi,M,L);
            TS_var = vc*p_var+vs*q_var;
            Psi = vc*p_Psi + vs*q_Psi;
            TS_Dvar = -vs*diag(1:M)*p_var(2:end,:) + vc(:,2:end)*diag(1:M)*q_var;

            res_sys = 0;
            res_var_phi = 0;
            % res_obs_phi = 0;
            res_target = 0;
            % res_PoincareCond = 0;
            for i = 1:N
                residue = FMAM_ODE.residue_system(system,Psi,TS_var,TS_Dvar,parameters,i);
                res_temp = fft(residue)/L;
                res_sys = max(res_sys,max(abs(res_temp(1:M+1))));

                [needMax, needMin] = obj.checkTarget(i);

                if needMax
                    p = p_var(:,i);
                    q = q_var(:,i);
                    phi = obj.stat.varPhiMax(i);

                    res_var_phi = max(res_var_phi,abs(FMAM_ODE.residue_phi_var(p,q,phi)));
                end

                if needMin
                    p = p_var(:,i);
                    q = q_var(:,i);
                    phi = obj.stat.varPhiMin(i);

                    res_var_phi = max(res_var_phi,abs(FMAM_ODE.residue_phi_var(p,q,phi)));
                end
            end

            observable = obj.obs;
            for j = 1:obj.n_per
                item_per = obj.items_perturb(j);
                propname = item_per.prop;

                idx_target = item_per.idx;
                if length(idx_target) > 1
                    idx_target = num2cell(idx_target);
                    idx_target = sub2ind(size(obj.stat.(propname)),idx_target{:});
                end
                val_target = obj.stat.(propname)(idx_target);
                if strcmpi(propname,'varAmp')
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi_max = obj.stat.varPhiMax(idx_target);
                    Phi_min = obj.stat.varPhiMin(idx_target);

                    res_temp = 2*val_target - (obj.var_target_curr(p,q,Phi_max) ...
                        - obj.var_target_curr(p,q,Phi_min));
                
                elseif strcmpi(propname,'varMax')
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi = obj.stat.varPhiMax(idx_target);

                    res_temp = val_target - obj.var_target_curr(p,q,Phi);
                
                elseif strcmpi(propname,'varMin')
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi = obj.stat.varPhiMin(idx_target);

                    res_temp = val_target - obj.var_target_curr(p,q,Phi);

                elseif strcmpi(propname,'obsAmp')
                    Phi_max = obj.stat.obsPhiMax(idx_target);
                    Phi_min = obj.stat.obsPhiMin(idx_target);

                    res_temp = 2*val_target - (obj.obs_target_curr(observable,p_var,q_var,Phi_max,idx_target) - ...
                        obj.obs_target_curr(observable,p_var,q_var,Phi_min,idx_target));

                elseif strcmpi(propname,'obsMax')
                    Phi = obj.stat.obsPhiMax(idx_target);
                    res_temp = val_target - obj.obs_target_curr(observable,p_var,q_var,Phi,idx_target);

                elseif strcmpi(propname,'obsMin')
                    Phi = obj.stat.obsPhiMin(idx_target);
                    res_temp = val_target - obj.obs_target_curr(observable,p_var,q_var,Phi,idx_target);

                elseif strcmpi(propname,'varPhase')
                    idx_target = item_per.idx;
                    idx_var1 = idx_target(1);
                    idx_var2 = idx_target(2);
                    idx_target = num2cell(idx_target);
                    idx_target = sub2ind(size(obj.stat.(propname)),idx_target{:});
                    val_target = obj.stat.(propname)(idx_target);

                    Phi1 = obj.stat.varPhiMax(idx_var1);
                    Phi2 = obj.stat.varPhiMax(idx_var2);

                    res_temp = val_target - FMAM_ODE.Trintegration(p_Psi,q_Psi,Phi1,Phi2);
                elseif strcmpi(propname,'obsPhase')

                else
                    res_temp = 0;
                end
                res_target = max(res_target,abs(res_temp));
            end

            val = [res_sys,res_var_phi,res_target];
        end

        function oneIter(obj)
            p_Psi_init_ = obj.p_Psi_init;
            q_Psi_init_ = obj.q_Psi_init;

            system = obj.sys;
            observable = obj.obs;
            Derivative_var = obj.derivatives.var;
            Derivative_obs = obj.derivatives.obs;
            DDerivative_obs = obj.derivatives.obs2;
            L = obj.Lconst;
            N = obj.dimVar;
            n = obj.dimObs;
            M = obj.truncationOrder;
            m = length(obj.stat.params);

            parameters = obj.stat.params;
            Coe_Controlled = obj.items_controlled;

            p_var = obj.stat.p_var;
            q_var = obj.stat.q_var;
            p_Psi = obj.stat.p_Psi;
            q_Psi = obj.stat.q_Psi;

            phi = (0:L-1)'*2*pi/L;
            [vc,vs] = obj.Vec_CS(phi,M,L);
            TS_var = vc*p_var+vs*q_var;
            Psi = vc*p_Psi + vs*q_Psi;
            TS_Dvar = -vs*diag(1:M)*p_var(2:end,:) + vc(:,2:end)*diag(1:M)*q_var;

            unknowns = {'params','p_Psi','q_Psi','p_var','q_var','varPhiMax','varPhiMin','obsPhiMax','obsPhiMin'};
            numUnknown = m + (N+1)*(2*M+1) + (N+n)*2;
            A = zeros(numUnknown,numUnknown);
            res = zeros(numUnknown,1);
            
            indexMap = struct;
            lastIdx = 0;

            for i = 1:size(unknowns,2)
                propname = unknowns{i};
                rows = size(obj.stat.(propname),1);
                columns = size(obj.stat.(propname),2);
                indexMap.(propname) = reshape(lastIdx + (1:rows*columns),[],columns);
                lastIdx = lastIdx + rows*columns;
            end

            % indexMap.params = 1:m;
            % lastIdx = m;
            % 
            % indexMap.tphi = lastIdx + 1;
            % lastIdx = indexMap.tphi;
            % 
            % indexMap.p_var = reshape(lastIdx + (1:N*(M+1)),[],N);
            % lastIdx = indexMap.p_var(end,end);
            % 
            % indexMap.q_var = reshape(lastIdx + (1:N*M),[],N);
            % lastIdx = indexMap.q_var(end,end);
            % 
            % indexMap.varPhiMax = lastIdx + (1:N);
            % lastIdx = indexMap.varPhiMax(end);
            % 
            % indexMap.varPhiMin = lastIdx + (1:N);
            % lastIdx = indexMap.varPhiMin(end);
            % 
            % indexMap.obsPhiMax = lastIdx + (1:n);
            % lastIdx = indexMap.obsPhiMax(end);
            % 
            % indexMap.obsPhiMin = lastIdx + (1:n);

            % Incremential equations of the ODE, num = N*(2*M+1)
             
            idxCollum = 1:m+(N+1)*(2*M+1);
            idxLastRow = 0;
            for i = 1:N
                [coe_params,coe_p_Psi,coe_q_Psi,coe_p_var,coe_q_var] = ...
                    FMAM_ODE.delta_coe_system(system,Derivative_var,Psi,TS_var, ...
                    parameters,vc,vs,i);
                A_temp = fft([coe_params,coe_p_Psi,coe_q_Psi,coe_p_var,coe_q_var])/L;
                A_plus1 = real(A_temp);
                A_plus2 = imag(A_temp);

                residue = FMAM_ODE.residue_system(system,Psi,TS_var,TS_Dvar,parameters,i);
                res_temp = fft(residue)/L;
                res_plus1 = real(res_temp);
                res_plus2 = imag(res_temp);

                A(idxLastRow + (1:M+1), idxCollum) = A_plus1(1:M+1,:);
                res(idxLastRow + (1:M+1)) = res_plus1(1:M+1);
                A(idxLastRow + M+1 + (1:M), idxCollum) = A_plus2(2:M+1,:);
                res(idxLastRow + M+1 + (1:M)) = res_plus2(2:M+1);

                idxLastRow = idxLastRow + (2*M+1);

                PV = obj.stat.PV;
 
                if obj.isPsiUpdated
                    if all([strcmpi(PV.name,'var'), (i == PV.idx)])
                        A(idxLastRow + (1:2), idxCollum) = A_plus1(M+2:M+3,:);
                        res(idxLastRow + (1:2)) = res_plus1(M+2:M+3);
                        idxLastRow = idxLastRow + 2;
                    elseif all([strcmpi(PV.name, 'obs'), (i == 1)])
                        A(idxLastRow + (1:2), idxCollum) = A_plus1(M+2:M+3,:);
                        res(idxLastRow + (1:2)) = res_plus1(M+2:M+3);
                        idxLastRow = idxLastRow + 2;
                    end
                else
                    A(idxLastRow + 1, idxCollum) = A_plus1(M+2,:);
                    res(idxLastRow + 1) = res_plus1(M+2);
                    idxLastRow = idxLastRow + 1;
                end
                
            end

            % num_extraParams = size(obj.items_controlled,2) - obj.n_per;
            % if obj.isPsiUpdated
            %     A(idxLastRow + (1:2+num_extraParams), idxCollum) = A_plus1(M+2:M+3+num_extraParams,:);
            %     res(idxLastRow + (1:2+num_extraParams)) = res_plus1(M+2:M+3+num_extraParams);
            %     idxLastRow = idxLastRow + 2 + num_extraParams;
            % else
            %     A(idxLastRow + (1:1+num_extraParams), idxCollum) = A_plus1(M+2:M+2+num_extraParams,:);
            %     res(idxLastRow + 1+num_extraParams) = res_plus1(M+2+num_extraParams);
            %     idxLastRow = idxLastRow + 1 + num_extraParams;
            % end

            % Incremental equations of extreme points, num = 2*(N+n)
            % Extreme points of variables

            for i = 1:N
                [needMax, needMin] = obj.checkTarget(i);
                idx_phi_max = indexMap.varPhiMax(i);
                idx_phi_min = indexMap.varPhiMin(i);
                if needMax
                    p = p_var(:,i);
                    q = q_var(:,i);
                    idx_p_var_i = indexMap.p_var(:,i);
                    idx_q_var_i = indexMap.q_var(:,i);
                    phi = obj.stat.varPhiMax(i);

                    [coe_p_var_i,coe_q_var_i,coe_phi] = FMAM_ODE.delta_coe_phi_var(p,q,phi);

                    A(idxLastRow+1,idx_p_var_i) = coe_p_var_i;
                    A(idxLastRow+1,idx_q_var_i) = coe_q_var_i;
                    A(idxLastRow+1,idx_phi_max) = coe_phi;
                    res(idxLastRow+1) = -FMAM_ODE.residue_phi_var(p,q,phi);
                else
                    A(idxLastRow+1,idx_phi_max) = 1;
                    res(idxLastRow+1) = 0;
                end
                idxLastRow = idxLastRow+1;

                if needMin
                    p = p_var(:,i);
                    q = q_var(:,i);
                    idx_p_var_i = indexMap.p_var(:,i);
                    idx_q_var_i = indexMap.q_var(:,i);
                    phi = obj.stat.varPhiMin(i);

                    [coe_p_var_i,coe_q_var_i,coe_phi] = FMAM_ODE.delta_coe_phi_var(p,q,phi);

                    A(idxLastRow+1,idx_p_var_i) = coe_p_var_i;
                    A(idxLastRow+1,idx_q_var_i) = coe_q_var_i;
                    A(idxLastRow+1,idx_phi_min) = coe_phi;
                    res(idxLastRow+1) = -FMAM_ODE.residue_phi_var(p,q,phi);
                else
                    A(idxLastRow+1,idx_phi_min) = 1;
                    res(idxLastRow+1) = 0;
                end
                idxLastRow = idxLastRow + 1;
            end
            
            % Extreme points of observables
            for k = 1:n
                [needMax, needMin] = obj.checkTarget(k);
                idx_phi_max = indexMap.obsPhiMax(k);
                idx_phi_min = indexMap.obsPhiMin(k);

                if needMax
                    phi = obj.stat.obsPhiMax(k);
                    [coe_p_var,coe_q_var,coe_phi] = ...
                        FMAM_ODE.delta_coe_obs_phi(Derivative_obs,DDerivative_obs,p_var,q_var,phi,k);

                    for i = 1:N
                        idx_p_var_i = indexMap.p_var(:,i);
                        idx_q_var_i = indexMap.q_var(:,i);
                        A(idxLastRow+1,idx_p_var_i) = coe_p_var(:,i)';
                        A(idxLastRow+1,idx_q_var_i) = coe_q_var(:,i)';
                    end
                    A(idxLastRow+1,idx_phi_max) = coe_phi;
                    res(idxLastRow+1) = -FMAM_ODE.residue_phi_obs(Derivative_obs,p_var,q_var,phi,k);
                else
                    A(idxLastRow+1,idx_phi_max) = 1;
                    res(idxLastRow+1) = 0;
                end
                idxLastRow = idxLastRow+1;

                if needMin
                    phi = obj.stat.obsPhiMin(k);
                    [coe_p_var,coe_q_var,coe_phi] = ...
                        FMAM_ODE.delta_coe_obs_phi(Derivative_obs,DDerivative_obs,p_var,q_var,phi,k);

                    for i = 1:N
                        idx_p_var_i = indexMap.p_var(:,i);
                        idx_q_var_i = indexMap.q_var(:,i);
                        A(idxLastRow+1,idx_p_var_i) = coe_p_var(:,i)';
                        A(idxLastRow+1,idx_q_var_i) = coe_q_var(:,i)';
                    end
                    A(idxLastRow+1,idx_phi_min) = coe_phi;
                    res(idxLastRow+1) = -FMAM_ODE.residue_phi_obs(Derivative_obs,p_var,q_var,phi,k);
                else
                    A(idxLastRow+1,idx_phi_min) = 1;
                    res(idxLastRow+1) = 0;
                end
                idxLastRow = idxLastRow+1;
            end

            % Restriction equations to fix the uncontrolled params, num =
            % m - n_controlled

            idx_params_fix = 1:m;
            idx_params_fix(Coe_Controlled) = [];
            for i = idx_params_fix
                A(idxLastRow+1,i) = FMAM_ODE.coeConst;
                res(idxLastRow+1) = 0;
                idxLastRow = idxLastRow + 1;
            end

            % Restriction equations to achieve modulation targets, num = n_per
            for j = 1:obj.n_per
                item_per = obj.items_perturb(j);
                propname = item_per.prop;
                if any([strcmpi(propname,'p_Psi'),strcmpi(propname,'q_Psi'),strcmpi(propname,'p_var'), ...
                        strcmpi(propname,'q_var'),strcmpi(propname,'varPhiMax'), ...
                        strcmpi(propname,'obsPhiMax'),strcmpi(propname,'varPhiMin'), ...
                        strcmpi(propname,'obsPhiMin')])
                    idx_target = item_per.idx;
                    if length(idx_target) > 1
                        idx_target = num2cell(idx_target);
                        idx_target = sub2ind(size(obj.stat.(propname)),idx_target{:});
                    end
                    A(idxLastRow+1, indexMap.(propname)(idx_target)) = FMAM_ODE.coeConst;
                    res(idxLastRow+1) = 0;
                elseif strcmpi(propname,'varAmp')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi_max = obj.stat.varPhiMax(idx_target);
                    Phi_min = obj.stat.varPhiMin(idx_target);

                    idx_p_var_i = indexMap.p_var(:,idx_target);
                    idx_q_var_i = indexMap.q_var(:,idx_target);
                    idx_phi_max = indexMap.varPhiMax(idx_target);
                    idx_phi_min = indexMap.varPhiMin(idx_target);

                    [coe_p_var_max,coe_q_var_max,coe_phi_var_max] = FMAM_ODE.delta_coe_var_target(p,q,Phi_max);
                    [coe_p_var_min,coe_q_var_min,coe_phi_var_min] = FMAM_ODE.delta_coe_var_target(p,q,Phi_min);

                    coe_p_var = coe_p_var_max - coe_p_var_min;
                    coe_q_var = coe_q_var_max - coe_q_var_min;
                    A(idxLastRow+1,idx_p_var_i) = coe_p_var;
                    A(idxLastRow+1,idx_q_var_i) = coe_q_var;
                    A(idxLastRow+1,idx_phi_max) = coe_phi_var_max;
                    A(idxLastRow+1,idx_phi_min) = -coe_phi_var_min;
                    res(idxLastRow+1) = 2*val_target - (FMAM_ODE.var_target_curr(p,q,Phi_max) - FMAM_ODE.var_target_curr(p,q,Phi_min));
                
                elseif strcmpi(propname,'varMax')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi = obj.stat.varPhiMax(idx_target);

                    idx_p_var_i = indexMap.p_var(:,idx_target);
                    idx_q_var_i = indexMap.q_var(:,idx_target);
                    idx_phi = indexMap.varPhiMax(idx_target);

                    [coe_p_var,coe_q_var,coe_phi_var] = FMAM_ODE.delta_coe_var_target(p,q,Phi);
                    A(idxLastRow+1,idx_p_var_i) = coe_p_var;
                    A(idxLastRow+1,idx_q_var_i) = coe_q_var;
                    A(idxLastRow+1,idx_phi) = coe_phi_var;
                    res(idxLastRow+1) = val_target - FMAM_ODE.var_target_curr(p,q,Phi);
                
                elseif strcmpi(propname,'varMin')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    p = p_var(:,idx_target);
                    q = q_var(:,idx_target);
                    Phi = obj.stat.varPhiMin(idx_target);

                    idx_p_var_i = indexMap.p_var(:,idx_target);
                    idx_q_var_i = indexMap.q_var(:,idx_target);
                    idx_phi = indexMap.varPhiMin(idx_target);

                    [coe_p_var,coe_q_var,coe_phi_var] = FMAM_ODE.delta_coe_var_target(p,q,Phi); 
                    A(idxLastRow+1,idx_p_var_i) = coe_p_var;
                    A(idxLastRow+1,idx_q_var_i) = coe_q_var;
                    A(idxLastRow+1,idx_phi) = coe_phi_var;
                    res(idxLastRow+1) = val_target - FMAM_ODE.var_target_curr(p,q,Phi);

                elseif strcmpi(propname,'obsAmp')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    Phi_max = obj.stat.obsPhiMax(idx_target);
                    Phi_min = obj.stat.obsPhiMin(idx_target);

                    idx_phi_max = indexMap.varPhiMax(idx_target);
                    idx_phi_min = indexMap.varPhiMin(idx_target);

                    [coe_p_var_max,coe_q_var_max,coe_phi_obs_max] = FMAM_ODE.delta_coe_obs_target ...
                        (Derivative_obs,p_var,q_var,Phi_max,idx_target);
                    [coe_p_var_min,coe_q_var_min,coe_phi_obs_min] = FMAM_ODE.delta_coe_obs_target ...
                        (Derivative_obs,p_var,q_var,Phi_min,idx_target);
                    for i = 1:N
                        idx_p_var_i = indexMap.p_var(:,i);
                        idx_q_var_i = indexMap.q_var(:,i);
                        A(idxLastRow+1,idx_p_var_i) = coe_p_var_max(:,i)' - coe_p_var_min(:,i)';
                        A(idxLastRow+1,idx_q_var_i) = coe_q_var_max(:,i)' - coe_q_var_min(:,i)';
                    end
                    A(idxLastRow+1,idx_phi_max) = coe_phi_obs_max;
                    A(idxLastRow+1,idx_phi_min) = -coe_phi_obs_min;
                    res(idxLastRow+1) = 2*val_target - (FMAM_ODE.obs_target_curr(observable,p_var,q_var,Phi_max,idx_target) - ...
                        FMAM_ODE.obs_target_curr(observable,p_var,q_var,Phi_min,idx_target));

                elseif strcmpi(propname,'obsMax')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    Phi = obj.stat.obsPhiMax(idx_target);
                    idx_phi = indexMap.obsPhiMax(idx_target);

                    [coe_p_var,coe_q_var,coe_phi_obs] = FMAM_ODE.delta_coe_obs_target ...
                        (Derivative_obs,p_var,q_var,Phi,idx_target);
                    for i = 1:N
                        idx_p_var_i = indexMap.p_var(:,i);
                        idx_q_var_i = indexMap.q_var(:,i);
                        A(idxLastRow+1,idx_p_var_i) = coe_p_var(:,i)';
                        A(idxLastRow+1,idx_q_var_i) = coe_q_var(:,i)';
                    end
                    A(idxLastRow+1,idx_phi) = coe_phi_obs;
                    res(idxLastRow+1) = val_target - FMAM_ODE.obs_target_curr(observable,p_var,q_var,Phi,idx_target);

                elseif strcmpi(propname,'obsMin')
                    idx_target = item_per.idx;
                    val_target = obj.stat.(propname)(idx_target);
                    Phi = obj.stat.obsPhiMin(idx_target);
                    idx_phi = indexMap.obsPhiMin(idx_target);

                    [coe_p_var,coe_q_var,coe_phi_obs] = FMAM_ODE.delta_coe_obs_target ...
                        (Derivative_obs,p_var,q_var,Phi,idx_target);
                    for i = 1:N
                        idx_p_var_i = indexMap.p_var(:,i);
                        idx_q_var_i = indexMap.q_var(:,i);
                        A(idxLastRow+1,idx_p_var_i) = coe_p_var(:,i)';
                        A(idxLastRow+1,idx_q_var_i) = coe_q_var(:,i)';
                    end
                    A(idxLastRow+1,idx_phi) = coe_phi_obs;
                    res(idxLastRow+1) = val_target - FMAM_ODE.obs_target_curr(observable,p_var,q_var,Phi,idx_target);

                elseif any([strcmpi(propname,'p_var_origin') strcmpi(propname,'q_var_origin')])

                elseif strcmpi(propname,'varPhase')
                    idx_target = item_per.idx;
                    idx_var1 = idx_target(1);
                    idx_var2 = idx_target(2);
                    
                    idx_target = num2cell(idx_target);
                    idx_target = sub2ind(size(obj.stat.(propname)),idx_target{:});
                    val_target = obj.stat.(propname)(idx_target);

                    Phi1 = obj.stat.varPhiMax(idx_var1);
                    Phi2 = obj.stat.varPhiMax(idx_var2);

                    idx_p_Psi = indexMap.p_Psi;
                    idx_q_Psi = indexMap.q_Psi;
                    idx_phi1 = indexMap.varPhiMax(idx_var1);
                    idx_phi2 = indexMap.varPhiMax(idx_var2);

                    [coe_p_Psi,coe_q_Psi,coe_phi1,coe_phi2] = FMAM_ODE.delta_coe_state_phase(p_Psi,q_Psi,Phi1,Phi2);

                    A(idxLastRow+1,idx_p_Psi) = coe_p_Psi;
                    A(idxLastRow+1,idx_q_Psi) = coe_q_Psi;
                    A(idxLastRow+1,idx_phi1) = coe_phi1;
                    A(idxLastRow+1,idx_phi2) = coe_phi2;
                    res(idxLastRow+1) = val_target - FMAM_ODE.Trintegration(p_Psi,q_Psi,Phi1,Phi2);
                elseif strcmpi(propname,'obsPhase')

                end
                idxLastRow = idxLastRow + 1;
            end
            % % Restriction equation to eliminate the degree of freedom in the
            % % tangent direction, num = 1
            % x0_curr = sum(p_var,1);
            % for i = 1:N
            %     fi = system{1,i};
            %     fi0_last = fi(x0_last,params_last);
            %     idx_p_var_i = indexMap.p_var(:,i);
            %     A(idxLastRow+1,idx_p_var_i) = fi0_last*ones(1,M+1);
            %     res(idxLastRow+1) = res(idxLastRow+1) - (x0_curr(i)-x0_last(i))*fi0_last;
            % end
            % idxLastRow = idxLastRow + 1;

            % Restriction equation for PV
            if obj.isPsiUpdated
                if strcmpi(PV.name,'var')
                    idx_PV = PV.idx;
                    idx_p = indexMap.p_var(3:end,idx_PV);
                    idx_q = indexMap.q_var(:,idx_PV);

                    A(idxLastRow + (1:M-1),idx_p) = FMAM_ODE.coeConst*diag(ones(M-1,1));
                    res(idxLastRow + (1:M-1)) = -p_var(3:end,idx_PV);

                    idxLastRow = idxLastRow + M-1;

                    A(idxLastRow + (1:M),idx_q) = FMAM_ODE.coeConst*diag(ones(M,1));
                    res(idxLastRow + (1:M)) = -q_var(:,idx_PV);
                elseif strcmpi(PV.name,'obs')
                    [coe_p_var,coe_q_var] = FMAM_ODE.delta_coe_observable(Derivative_obs,p_var,q_var,PV);
                    idx_column = m+2*M+1 + (1:N*(2*M+1));
                    A_temp = fft([coe_p_var,coe_q_var])/L;
                    A_plus1=real(A_temp);
                    A_plus2=imag(A_temp);
    
                    f = observable{PV.idx};
                    res_obs = f(TS_var);
                    res_temp = fft(res_obs)/L;
                    res_plus1 = real(res_temp);
                    res_plus2 = imag(res_temp);
    
                    A(idxLastRow + (1:M-1),idx_column) = A_plus1(3:M+1,:);
                    res(idxLastRow + (1:M-1)) = -res_plus1(3:M+1);
                    idxLastRow = idxLastRow + M-1;
    
                    A(idxLastRow + (1:M),idx_column) = A_plus2(2:M+1,:);
                    res(idxLastRow + (1:M)) = -res_plus2(2:M+1);
                end
            else
                idx_p_Psi = indexMap.p_Psi(2:end);
                idx_q_Psi = indexMap.q_Psi;
                A(idxLastRow + (1:M),idx_p_Psi) = diag(ones(M,1));
                res(idxLastRow + (1:M)) = -(p_Psi(2:end) - p_Psi_init_(2:end));

                idxLastRow = idxLastRow + M;

                A(idxLastRow + (1:M),idx_q_Psi) = diag(ones(M,1));
                res(idxLastRow + (1:M)) = -(q_Psi - q_Psi_init_);
            end
            
            % Solve the system of linear incremential equations
            increments = A\res;

            % Update obj.stat
            for i = 1:size(unknowns,2)
                propname = unknowns{i};
                idxProp = indexMap.(propname);
                columns = size(idxProp,2);
                propIncrement = reshape(increments(idxProp(:)),[],columns);
                obj.stat.add(propname,'all',obj.lr*propIncrement)
            end
            if max(obj.res()) > obj.maxerr
               for i = 1:size(unknowns,2)
                   propname = unknowns{i};
                   idxProp = indexMap.(propname);
                   columns = size(idxProp,2);
                   propIncrement = reshape(increments(idxProp(:)),[],columns);
                   obj.stat.minus(propname,'all',obj.lr*propIncrement)
               end
               obj.lr = max(obj.lr/2,obj.minLr);
            else
               obj.lr = min(obj.lr*2,obj.maxLr);
            end

            % Update obj.stat

            % for i = 1:size(unknowns,2)
            %     propname = unknowns{i};
            %     idxProp = indexMap.(propname);
            %     columns = size(idxProp,2);
            %     propIncrement = reshape(increments(idxProp(:)),[],columns);
            %     obj.stat.add(propname,'all',FMAM_ODE.lr*propIncrement)
            % end

            if max(obj.res()) > obj.maxerr
                disp(increments)
            end
        end    

        function [needMax,needMin] = checkTarget(obj,idx)
            needMax = false;
            needMin = false;
            items = obj.items_perturb;
            numPerturb = size(items,2);

            for j = 1:numPerturb
                if items(j).idx == idx
                    propname = items(j).prop;
                    needMax = any([...
                        strcmpi(propname,'varMax'),...
                        strcmpi(propname,'varPhase'),...
                        strcmpi(propname,'varAmp')
                        ]);
                    needMin = any([...
                        strcmpi(propname,'varMin'),...
                        strcmpi(propname,'varAmp')
                        ]);
                end
                
                if all([needMax,needMin])
                    break
                end
            end
        end

    end

    methods(Static)
        function [vc,vs]=Vec_CS(phi,M,L)
            vc=zeros(L,M+1);
            vs=zeros(L,M);
            for i= 1:M
                vc(:,i)=cos((i-1)*phi);
                vs(:,i)=sin(i*phi);
            end
            vc(:,M+1)=cos(M*phi);
        end

        function output = Trintegration(p,q,phi_L,phi_U)
            % Compute the integration of a trigonometric polynomial with coefficients
            % saved in p and q. phi_L and phi_U are the lower bound and upper bound of
            % the integral, respectively.
            M = size(q,1);
            [vc_U,vs_U] = FMAM_ODE.Vec_CS(phi_U,M,1);
            [vc_L,vs_L] = FMAM_ODE.Vec_CS(phi_L,M,1);
            U = p(1)*phi_U+vs_U*diag((1:M).^(-1))*p(2:end)-vc_U(:,2:end)*diag((1:M).^(-1))*q;
            L = p(1)*phi_L+vs_L*diag((1:M).^(-1))*p(2:end)-vc_L(:,2:end)*diag((1:M).^(-1))*q;
            output = U-L;
        end

        function [coe_para,coe_p_Psi,coe_q_Psi,coe_p_var,coe_q_var] = ...
         delta_coe_system(System,Derivative,tphi,TS_variable,parameters,...
         vc,vs,idx_eq)
            M = size(vs,2);
            L = size(TS_variable,1);
            N = size(TS_variable,2);
            m = size(parameters,2);

            coe_para = zeros(L,m);

            F = System{1,idx_eq};
        
            for i = 1:m
                coe_para(:,i) = -tphi.*Derivative(idx_eq,N+i).function(TS_variable,parameters);
            end
            coe_p_Psi = -diag(F(TS_variable,parameters))*vc;       
            coe_q_Psi = -diag(F(TS_variable,parameters))*vs; 

            A = zeros(L,N);
            for i = 1:N 
                A(:,i) = -tphi.*Derivative(idx_eq,i).function(TS_variable,parameters);
            end

            coe_p_var = zeros(L,N*(M+1));
            coe_q_var = zeros(L,N*M);
            for i = 1:N
                idx_p = (i-1)*(M+1) + (1:M+1);
                idx_q = (i-1)*M + (1:M);
                coe_p_var(:,idx_p) = diag(A(:,i))*vc-(idx_eq==i)*[zeros(size(vs,1),1),vs*diag(1:M)];
                coe_q_var(:,idx_q) = diag(A(:,i))*vs+(idx_eq==i)*vc(:,2:end)*diag(1:M);
            end
        end

        function output = residue_system(System,tphi,TS_variable, ...
            TS_Dvariable,parameters,idx_eq)
            F = System{idx_eq};
            output = -TS_Dvariable(:,idx_eq) + tphi.*F(TS_variable,parameters);
        end

        function [coe_p_var,coe_q_var,coe_phi] = delta_coe_phi_var(p,q,phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
            coe_p_var = -(0:M).*[0,vs];
            coe_q_var = (1:M).*vc(2:end);
            coe_phi = -((0:M).^2.*vc*p+(1:M).^2.*vs*q);
        end

        function [coe_p_var,coe_q_var,coe_phi] = delta_coe_obs_phi(Derivative_obs, ...
                DDerivative_obs, p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            coe_phi = 0;
            coe_p_var = zeros(size(p_var));
            coe_q_var = zeros(size(q_var));
        
            pt_var = vc*p_var+vs*q_var; % evaluated at a single time point
            pt_Dvar = -vs*diag(1:M)*p_var(2:end,:) + vc(:,2:end)*diag(1:M)*q_var;
            for i = 1:N
                dxidphi = pt_Dvar(i);
                dykdxi = Derivative_obs(k,i).function(pt_var);
        
                coe_p_var(:,i) = coe_p_var(:,i) + (-dykdxi*(0:M).*[0 vs])';
                coe_q_var(:,i) = coe_q_var(:,i) + (dykdxi*(1:M).*vc(2:end))';
                coe_phi = coe_phi-dykdxi*((0:M).^2.*vc*p_var(:,i)+(1:M).^2.*vs*q_var(:,i));
                for j = 1:N
                    dxjdphi = pt_Dvar(j);
                    ddykdxixj = DDerivative_obs(k,i,j).function(pt_var);
                    coe_p_var(:,j) = coe_p_var(:,j) + ddykdxixj*dxidphi*vc';
                    coe_q_var(:,j) = coe_q_var(:,j) + ddykdxixj*dxidphi*vs';
                    coe_phi = coe_phi + dxidphi*ddykdxixj*dxjdphi;
                end
            end
        end

        function res = residue_phi_var(p,q,phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
            res = -vs*diag(1:M)*p(2:end,:) + vc(:,2:end)*diag(1:M)*q;
        end

        function res = residue_phi_obs(Derivative_obs,p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            res = 0;
            pt_var = vc*p_var+vs*q_var; % evaluated at a single time point
            pt_Dvar = -vs*diag(1:M)*p_var(2:end,:) + vc(:,2:end)*diag(1:M)*q_var;
            for i = 1:N
                dykdxi = Derivative_obs(k,i).function(pt_var);
                dxidphi = pt_Dvar(i);
                res = res-dykdxi*dxidphi;
            end
        end

        function [coe_p_var,coe_q_var,coe_phi_var] = delta_coe_var_target(p,q,Phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);
        
            coe_p_var = vc;
            coe_q_var = vs;
        
            coe_phi_var = (1:M).*vc(2:end)*q-(0:M).*[0 vs]*p;
        end

        function [coe_p_var,coe_q_var] = delta_coe_observable(Derivative_obs,p_var,q_var,PV)
            assert(strcmpi(PV.name,'obs'))
            N = size(q_var,2);
            M = size(q_var,1);
            L = FMAM_ODE.Lconst;
            phi = (0:L-1)'*2*pi/L;

            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,L);

            TS_var = vc*p_var + vs*q_var;
        
            A = zeros(L,N);
            for i = 1:N
                temp = -Derivative_obs(PV.idx,i).function(TS_var);
                A(:,i) = ones(L,1).*temp;
            end

            coe_p_var = zeros(L,N*(M+1));
            coe_q_var = zeros(L,N*M);
            for i = 1:N
                idx_p = (i-1)*(M+1) + (1:M+1);
                idx_q = (i-1)*M + (1:M);
                coe_p_var(:,idx_p) = diag(A(:,i))*vc;
                coe_q_var(:,idx_q) = diag(A(:,i))*vs;
            end
        end

        function output = var_target_curr(p,q,Phi)
            M = size(q,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);
            
            output = vc*p+vs*q;
        end

        function [coe_p_var,coe_q_var,coe_phi_obs] = ...
                 delta_coe_obs_target(Derivative_obs,p_var,q_var,phi,k)
            M = size(q_var,1);
            N = size(q_var,2);
            [vc,vs] = FMAM_ODE.Vec_CS(phi,M,1);
        
            coe_p_var = zeros(size(p_var));
            coe_q_var = zeros(size(q_var));
            coe_phi_obs = 0;
        
            pt_var = vc*p_var+vs*q_var;
            pt_Dvar = -vs*diag(1:M)*p_var(2:end,:) + vc(:,2:end)*diag(1:M)*q_var;
            for i = 1:N
                dykdxi = Derivative_obs(k,i).function(pt_var);
                dxidphi = pt_Dvar(i);

                coe_p_var(:,i) = dykdxi*vc';
                coe_q_var(:,i) = dykdxi*vs';
                coe_phi_obs = coe_phi_obs + dykdxi*dxidphi;
            end
        end

        function output = obs_target_curr(obs,p_var,q_var,Phi,k)
            M = size(q_var,1);
            [vc,vs] = FMAM_ODE.Vec_CS(Phi,M,1);  
            pt_var = vc*p_var+vs*q_var;
        
            F = obs{k};
            output = F(pt_var);
        end

        function [coe_p_Psi,coe_q_Psi,coe_phi1,coe_phi2] = delta_coe_state_phase(p,q,Phi1,Phi2)
            M = size(q,1);
            [vc1,vs1] = FMAM_ODE.Vec_CS(Phi1,M,1);
            [vc2,vs2] = FMAM_ODE.Vec_CS(Phi2,M,1);
            coe_p_Psi = [Phi2-Phi1 (vs2-vs1)*diag((1:M).^(-1))];
            coe_q_Psi = -(vc2(:,2:end)-vc1(:,2:end))*diag((1:M).^(-1));
            coe_phi1 = -(vc1*p+vs1*q);
            coe_phi2 = vc2*p+vs2*q;
        end
    end
end

