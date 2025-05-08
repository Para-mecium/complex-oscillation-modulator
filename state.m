classdef state < handle
    properties (Constant)
        Lconst = 20000
        LphiConst = 500
        countMax = 20
        errMax = 1e-2
    end

    properties
        PV
    end

    properties (SetAccess=private)
        obs
        % category 1
        params
        p_Psi
        q_Psi
        p_var
        q_var
        p_var_origin
        q_var_origin
        p_obs
        q_obs
        % category 2
        a
        b
        varPhiMax = []
        obsPhiMax = []
        varPhiMin = []
        obsPhiMin = []
        varAmp = []
        obsAmp = []
        period = []
        varMax = []
        obsMax = []
        varMin = []
        obsMin = []
        varPhase = []
        obsPhase = []
    end

    properties(Dependent)
        dimSys
        dimObs
        dimParams
        truncationOrder
        Atrans

        t
        phi
        Psi
        TS_var
        TS_obs
    end
    methods
        function obj = state(obs,Params,t,TS_var,M,PV)
            obj.obs = obs;
            obj.params = Params;
            TS_obs = state.getObs(obs,TS_var);
            [a_PV,b_PV,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = state.fourierCoeffs(M,t,TS_var,TS_obs,PV);
            obj.a = a_PV;
            obj.b = b_PV;
            obj.p_Psi = p_Psi;
            obj.q_Psi = q_Psi;
            obj.p_var = p_variable;
            obj.q_var = q_variable;
            obj.p_obs = p_observable;
            obj.q_obs = q_observable;
            obj.period = 2*pi*p_Psi(1);
            obj.PV = PV;

            obj.updateFCOrigin();
            obj.updateVar2();
            obj.updateObs2();
        end
        
        function prop = get.dimSys(obj)
            prop = size(obj.p_var,2);
        end

        function prop = get.dimObs(obj)
            prop = size(obj.p_obs,2);
        end

        function prop = get.dimParams(obj)
            prop = size(obj.params,2);
        end

        function prop = get.truncationOrder(obj)
            prop = size(obj.q_var,1);
        end

        function prop = get.TS_var(obj)
            L = obj.LphiConst;
            M = obj.truncationOrder;

            [vc,vs] = FMAM_ODE.Vec_CS(obj.phi,M,L);
            prop = vc*obj.p_var+vs*obj.q_var;
        end

        function prop = get.TS_obs(obj)
            L = obj.LphiConst;

            funcs = obj.obs;
            prop = zeros(L,obj.dimObs);
            for k = 1:obj.dimObs
                func = funcs{k};
                prop(:,k) = func(obj.TS_var);
            end
        end

        function prop = get.Psi(obj)
            L = obj.LphiConst;
            M = obj.truncationOrder;

            [vc,vs] = FMAM_ODE.Vec_CS(obj.phi,M,L);
            prop = vc*obj.p_Psi+vs*obj.q_Psi;
        end

        function prop = get.phi(obj)
            L = obj.LphiConst;
            prop = (0:L-1)'*2*pi/L;
        end

        function prop = get.Atrans(obj)
            L = obj.LphiConst;
            prop = zeros(L,L+1);
            prop(1,1) = -3;prop(1,2) = 4;prop(1,3) = -1;
            for j = 2:L
                prop(j,j+1) = 1;
                prop(j,j-1) = -1;
            end
            prop = prop/(4*pi/L);
        end

        function prop = get.t(obj)
            Atrans_ = obj.Atrans;
            Atrans1 = Atrans_(:,2:end);

            t_increment = Atrans1\obj.Psi;
            prop = [0;t_increment];
            prop = prop(1:end-1);
        end

        function updateVar2(obj)
            % update the variable properties of category 2 
            Phi_max = zeros(1,obj.dimSys);
            Phi_min = zeros(1,obj.dimSys);
            Amp = zeros(1,obj.dimSys);
            VariableMax = zeros(1,obj.dimSys);
            VariableMin = zeros(1,obj.dimSys);

            for i = 1:obj.dimSys
                [Phi_max(i),Phi_min(i),Amp(i),VariableMax(i),VariableMin(i)] = obj.FindExtreme(obj.p_var(:,i),obj.q_var(:,i),obj.LphiConst);
            end
            tMax = state.Trintegration(obj.p_Psi,obj.q_Psi,zeros(obj.dimSys,1),Phi_max');
            Phase = repmat(tMax',obj.dimSys,1) - repmat(tMax,1,obj.dimSys);

            obj.varPhiMax = Phi_max;
            obj.varPhiMin = Phi_min;
            obj.varAmp = Amp;
            obj.varMax = VariableMax;
            obj.varMin = VariableMin;
            obj.varPhase = Phase;
        end

        function updateObs2(obj)
            if obj.dimObs == 0
                return
            end
            % update the observable properties of category 2 
            Phi_max = zeros(1,obj.dimObs);
            Phi_min = zeros(1,obj.dimObs);
            amp = zeros(1,obj.dimObs);
            VariableMax = zeros(1,obj.dimObs);
            VariableMin = zeros(1,obj.dimObs);

            for i = 1:obj.dimObs
                [Phi_max(i),Phi_min(i),amp(i),VariableMax(i),VariableMin(i)] = obj.FindExtreme(obj.p_obs(:,i),obj.q_obs(:,i),obj.Lconst);
            end
            tMax = state.Trintegration(obj.p_Psi,obj.q_Psi,zeros(obj.dimObs,1),Phi_max');
            Phase = repmat(tMax',obj.dimObs,1) - repmat(tMax,1,obj.dimObs);

            obj.obsPhiMax = Phi_max;
            obj.obsPhiMin = Phi_min;
            obj.obsAmp = amp;
            obj.obsMax = VariableMax;
            obj.obsMin = VariableMin;
            obj.obsPhase = Phase;
        end

        function updatePeriod(obj)
            obj.period = 2*pi*obj.p_Psi(1);
        end

        function add(obj,propname,idx,val)
            if strcmpi(idx,'all')
                if ~isequal(size(obj.(propname)),size(val))
                    val = val';
                end
                obj.(propname) = obj.(propname) + val;
            else
                obj.(propname)(idx) = obj.(propname)(idx) + val;
            end
        end

        function minus(obj,propname,idx,val)
            if strcmpi(idx,'all')
                if ~isequal(size(obj.(propname)),size(val))
                    val = val';
                end
                obj.(propname) = obj.(propname) - val;
            else
                obj.(propname)(idx) = obj.(propname)(idx) - val;
            end
        end

        function changeTruncOrder(obj,M)
            M1 = obj.truncationOrder;
            if M1 <= M
                obj.p_var = obj.p_var(1:M1+1,:);
                obj.q_var = obj.q_var(1:M1,:);
            else
                obj.p_var = [obj.p_var;zeros(M-M1,obj.dimSys)];
                obj.q_var = [obj.q_var;zeros(M-M1,obj.dimSys)];
            end
        end

        function TSplot(obj,propname,idx,timedomain)
            figure
            grid on
            plot(obj.(timedomain),[obj.(propname)(:,idx)])

            xlabel(timedomain)
            ylabel([propname ', idx = ' num2str(idx)])
        end

        function updatePV(obj,PV)
            obj.PV = PV;
            M = obj.truncationOrder;
            t_ = obj.t;
            TS_var_ = obj.TS_var;
            TS_obs_ = obj.getObs(obj.obs,TS_var_);

            [p_Psi_,q_Psi_,p_var_,q_var_,p_obs_,q_obs_] = ...
            state.fourierCoeffs(M,t_,TS_var_,TS_obs_,PV);

            obj.p_Psi = p_Psi_;
            obj.q_Psi = q_Psi_;
            obj.p_var = p_var_;
            obj.q_var = q_var_;
            obj.p_obs = p_obs_;
            obj.q_obs = q_obs_;

            obj.updateVar2();
            obj.updateObs2();
            obj.updatePeriod();
        end

        function updateFCOrigin(obj)
            N = obj.dimSys;
            L = obj.LphiConst;
            M = obj.truncationOrder;
            T = obj.p_Psi(1)*2*pi;
            
            t_ = obj.t/(T/(2*pi));
            TS_var_ = obj.TS_var;
            t_equal = (0:L-1)'/L*2*pi;

            for i=1:N
                TS_var_(:,i) = spline(t_,TS_var_(:,i),t_equal);
            end
            temp = fft(TS_var_);
            p_variable = real(temp)/L;
            q_variable = imag(temp)/L;

            p_variable = 2*p_variable(1:M+1,:);
            p_variable(1,:) = p_variable(1,:)/2;
            q_variable = -2*q_variable(2:M+1,:);

            obj.p_var_origin = p_variable;
            obj.q_var_origin = q_variable;
        end
    end

    
    methods(Static)
        function TS_obs = getObs(obs,TS_var)
            TS_obs = zeros(size(TS_var,1),size(obs,2));
            for k = 1:size(obs,2)
                funcObs = obs{1,k};
                TS_obs(:,k) = funcObs(TS_var);
            end
        end
        function [phi_max,phi_min,amplitude,var_max,var_min] = FindExtreme(p,q,L)
            M = size(q,1);
            count = 0;
            while count < state.countMax
                phi = (0:L-1)'/L*2*pi;
                [vc,vs] = FMAM_ODE.Vec_CS(phi,M,L);
                
                TS = vc*p+vs*q;

                [var_max,index_max] = max(TS);
                [var_min,index_min] = min(TS);

                phi_max = phi(index_max);
                phi_min = phi(index_min);
    
                err_phiMax = abs(FMAM_ODE.residue_phi_var(p,q,phi_max));
                err_phiMin = abs(FMAM_ODE.residue_phi_var(p,q,phi_min));
                err = max(err_phiMax,err_phiMin);
                if err < state.errMax
                    break
                else 
                    L = L*2;
                end
            end

            phi_max = phi(index_max);
            phi_min = phi(index_min);
            amplitude = (var_max-var_min)/2;
        end

        function output = Trintegration(p,q,phi_L,phi_U)
        % Compute the integration of a trigonometric polynomial with coefficients
        % saved in p and q. phi_L and phi_U are the lower bound and upper bound of
        % the integral, respectively.
            M = size(q,1);
            assert(isequal(size(phi_L),size(phi_U)))
            ptNum = length(phi_L);
            [vc_U,vs_U] = FMAM_ODE.Vec_CS(phi_U,M,ptNum);
            [vc_L,vs_L] = FMAM_ODE.Vec_CS(phi_L,M,ptNum);
            U = p(1)*phi_U+vs_U*diag((1:M).^(-1))*p(2:end)-vc_U(:,2:end)*diag((1:M).^(-1))*q;
            L = p(1)*phi_L+vs_L*diag((1:M).^(-1))*p(2:end)-vc_L(:,2:end)*diag((1:M).^(-1))*q;
            output = U-L;
        end

        function [a,b,p_Psi,q_Psi,p_variable,q_variable,p_observable,q_observable] = ...
            fourierCoeffs(M,t,TS_variable,TS_observable,PV)
            dim = size(TS_variable,2);
            n = size(TS_observable,2);
            p_observable = [];
            q_observable = [];
            
            T = t(end);

            L = state.Lconst;
            Lphi = state.LphiConst;
            phi = (0:Lphi-1)'*2*pi/Lphi;
            t1 = linspace(0,T,L)';
            TS_variable_1 = zeros(L,dim);
            TS_observable_1 = zeros(L,n);
            for i=1:dim
                TS_variable_1(:,i) = spline(t,TS_variable(:,i),t1);
            end
            
            for i = 1:n
                TS_observable_1(:,i) = spline(t,TS_observable(:,i),t1);
            end

            if strcmpi(PV.name,'var')
                X = TS_variable_1(:,PV.idx);
            elseif strcmpi(PV.name,'obs')
                X = TS_observable_1(:,PV.idx);
            else
                error('Please check the class of the parimary variable.')
            end
            
            [~,I1]=max(X);
            
            % Translation of the time in order that the regualtory variable starts at
            % the maximum value.
            TS_variable_1 = [TS_variable_1(I1:end-1,:);TS_variable_1(1:I1,:)];
            X = [X(I1:end-1);X(1:I1)];
            if n>0
                TS_observable_1 = [TS_observable_1(I1:end-1,:);TS_observable_1(1:I1,:)];
            end
            
            t1 = [t1(I1:end-1)-t1(I1);t1(1:I1)+T-t1(I1)];
            
            [~,I2] = min(X);
            a=(max(X)-min(X))/2;% a is a scalar, represents the initial amplitude.
            b=(max(X)+min(X))/2;% b is a scalar.

            %% Computation of the time rescaling function Psi
            % tt is the vector of the same length with phi. tt record the time
            % at which the regulatory variable equals to the corresponding value in the
            % phi-field. tt gives the functional relation between phi and t. Notice
            % that the curve of acos+b has two intersections with every horizontal line
            % between b-a and b+a. Therefore, [0,2*pi] should be separated into two
            % halves.

            % Difference matrix
            Atrans = zeros(Lphi,Lphi+1);
            Atrans(1,1) = -3;Atrans(1,2) = 4;Atrans(1,3) = -1;
            for j = 2:Lphi
                Atrans(j,j+1) = 1;
                Atrans(j,j-1) = -1;
            end
            Atrans = Atrans/(4*pi/Lphi);

            t_phi = zeros(Lphi,1);
            index = zeros(Lphi,1);
                for i = 1:Lphi
                    if i<=round(Lphi/2)
                        [~,I] = min(abs(a*cos(phi(i))+b-X(1:I2-1)));
                        t_phi(i) = t1(I);
                        index(i) = I;
                    else
                        [~,I] = min(abs(a*cos(phi(i))+b-X(I2+1:end)));
                        t_phi(i) = t1(I2+I);
                        index(i) = I2+I;
                    end
                end
            
            t = [t_phi;T];
            TS_variable = TS_variable_1(index,:);
            if n>0
                TS_observable = TS_observable_1(index,:);
            end

            Psi=Atrans*t;% compute the derivative of Psi to phi using difference method

            
            %% Transform the time sequences into Fourier series
            temp_Psi = fft(Psi);
            temp_variable = fft(TS_variable);
            
            if n > 0
                temp_observable = fft(TS_observable);
            end
            % The first collum of p saves the cosine coefficients of Psi, the rest collums 
            % save the cosine coefficients of the other variables.
            % The first collum of q saves the sine coefficients of Psi, the rest
            % collums save the sine coefficients of the other variables. It should be
            % noticed that p has one more row than q because of the constant term.
            p_Psi = real(temp_Psi)/Lphi;
            q_Psi = imag(temp_Psi)/Lphi;
            p_variable = real(temp_variable)/Lphi;
            q_variable = imag(temp_variable)/Lphi;
            if n>0
                p_observable = real(temp_observable)/Lphi;
                q_observable = imag(temp_observable)/Lphi;
            end
            
            % Truncate the Fourier series up to the frequency M/(2*pi).
            p_Psi=2*p_Psi(1:M+1,:);
            p_Psi(1,:)=p_Psi(1,:)/2;
            q_Psi=-2*q_Psi(2:M+1,:);
            p_variable=2*p_variable(1:M+1,:);
            p_variable(1,:)=p_variable(1,:)/2;
            q_variable=-2*q_variable(2:M+1,:);
            if n>0
                p_observable=2*p_observable(1:M+1,:);
                p_observable(1,:)=p_observable(1,:)/2;
                q_observable=-2*q_observable(2:M+1,:);
            end
            if strcmpi(PV.name,'var')
                p_variable(1,PV.idx) = b;
                p_variable(2,PV.idx) = a;
                p_variable(3:end,PV.idx) = zeros(M-1,1);
                q_variable(:,PV.idx) = zeros(M,1);
            elseif strcmpi(PV.name,'obs')
                p_observable(3:end,PV.idx) = zeros(M-1,1);
                q_observable(:,PV.idx) = zeros(M,1);
            end
        end
    end
end

