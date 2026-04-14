classdef SimConf_LyapExp < CompConf

    properties
        methodname
        defaultrefine
        priority
        allowEvents
        simconf;
        type;
        name;
        pointlabel = 'P';
        nrActiveMax = 1;
    end
   
    methods
        
        function s = getLabel(obj)
            s = ['point_' obj.type '_' obj.pointlabel '_' obj.simconf.getLabel()];
        end
        function b = hasTarget(~) 
            b = [];
        end
        function p = getPrioritynumber(obj)
            % p = getPrioritynumber@CompConf();
            p=-1;
        end
        function l = getSolutionLabel(obj)
            l = obj.type;%'LE';
        end
        function s = getName(obj)
            s = obj.name;%'LyapExp';
        end
        function s = toString(obj)
            s = sprintf('%s', obj.getName());
        end
        function s = toStringCL(obj) 
            s = sprintf('%s (%s)',obj.getName(), obj.methodname);
        end

        function oi = getOutputInterpreter(~)
            oi = CLSimOutputInterpreter();
        end

        function m = getMethodName(obj)
           m = obj.simconf.getMethodName(); 
        end
        
        function obj = SimConf_LyapExp(name, type, simconf)
            % obj = obj@SimConf_ConnectBase(name, type, simconf, @GUIConnect_getPointHet);  
            obj.simconf = simconf;
            obj.type = type;
            obj.name = name;
            obj.defaultrefine = simconf.defaultrefine;
            obj.methodname = simconf.methodname;
            obj.priority = simconf.priority;
            obj.allowEvents = simconf.allowEvents; 

           end


        function b = isAvailable(obj, settings)
            % change this to the default value to prevent the starter
            % window from being the first thing that opens when you open
            % the program. 
            initialpoint = settings.IP;
            b = ~isempty(settings.system) && strcmp(initialpoint.getILabel(), 'P');
            % b = 1;
            
        end

        function b = isHidden(obj)
            % change this to the default value to prevent the starter
            % window from being the first thing that opens when you open
            % the program. 
            % b = obj.simconf.isHidden();
            b = 1;
            
        end


        function configureSettings(obj, settings)
            configureSettings@CompConf(obj, settings);
            
            obj.install(settings, 'forward', true, InputRestrictions.BOOL, [0, 1, 1]);
            settings.setVisible('forward', false);
            %odemethods = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'ode78', 'ode87'};
            %obj.install(settings, 'TMP_Method_TMP', obj.methodname, InputRestrictionCategory(odemethods), [1, 1, 1]);
            
            system = settings.system;
        
            if ~isempty(system) % this means that we are dealing with an empty settings object, independent of a system 
                settings.addSetting('time', CLSetting(system.getTimeName(), 0, InputRestrictions.NUM, 2, 1, 0, CLSettingsHelp.getHelp('time') ));
                coordinates = CLSettingCoordinates(settings, system.getCoordinates());
                settings.addSetting('coord', coordinates);
        
                parameters = settings.getSetting('parameters');
                
                if isempty(parameters)
                    parameters = CLSettingParameters(settings, system.getParameters(),true);
                    settings.addSetting('parameters', parameters);
                else
                    parameters.revive(settings, true); %active: allow for one parameter to be an array
                end
                
                obj.install(settings, 'Interval', 1, InputRestrictions.NUM, [3, 1, 1]);
                
                settings.addSetting('InitStepSize_sim', CLSettingBlank('InitStepSize', [], InputRestrictions.POS, 3, 1, 3, CLSettingsHelp.getHelp('InitStepSize_sim'), '<automatic>'));
                settings.addSetting('MaxStepSize_sim', CLSettingBlank('MaxStepSize', [], InputRestrictions.POS, 3, 1, 4, CLSettingsHelp.getHelp('MaxStepSize_sim'), '<automatic>'));
                settings.addSetting('con_ortho',    CLSetting('Orthogonalization Interval',    1   , InputRestrictions.POS, 2, 8, 1, ''));
                settings.addSetting('con_storage',  CLSetting('Storage steps',  20   , InputRestrictions.INT_ge0, 2, 8, 2, ''));
                settings.addSetting('ReUsePrevComp',CLSetting('Restart from last state', 1, InputRestrictions.BOOL, 2, 8, 3, CLSettingsHelp.getHelp('ReUsePrevComp')));
                settings.addSetting('ParamVals',   CLSetting('Active Parameter Range',[]  , InputRestrictions.VECTOR, 2, 8, 4, CLSettingsHelp.getHelp('ReUsePrevComp')));

                obj.install(settings, 'RelTolerance', 1e-3, InputRestrictions.POS, [3, 1, 5]);
                obj.install(settings, 'AbsTolerance', 1e-6, InputRestrictions.POS, [3, 1, 6]);  %can also be vector (component wise) FIX?
                obj.install(settings, 'Refine', obj.defaultrefine, InputRestrictions.INT_g0, [3, 1, 7]);  %remark: ODE45 is '4' FIX?
                obj.install(settings, 'Normcontrol', false, InputRestrictions.BOOL, [3, 1, 8]);
                
                % We do not allow for Event functions, this will interfere
                % with speed too much.
                % settings.addSetting('eventfunction' , CLSettingEventFunction('EventFunction'));
                % ef = settings.getSetting('eventfunction');
                % ef.setVisible(obj.allowEvents);
                
                refinesetting = settings.getSetting('Refine');
                %if ~refinesetting.isAdjusted()
                refinesetting.forceValue(obj.defaultrefine);
            end
            
        end
        
        function alist = actions(obj)
            alist(1).label = 'Forward';
            alist(1).function = @(settings, solution) obj.compute(settings, true);
            alist(1).valid = @(~,~) 1;
            % We will not support Lyapunov Exponents in Backward time.
            % alist(2).label = 'Backward'; % TODO: implement the lyapunov
            % exponent calculator to allow for extending the curve
            % alist(2).function = @(settings, solution) obj.compute(settings, false);
            % alist(2).valid = @(~,~) 1;
            % 
            % alist(3).label = '-'; alist(3).function = [];
            % 
            
        end
        
        function [solution, errmsg, overwrite] = compute(obj, settings, forward)
            solution = []; errmsg = []; overwrite = 0;
            global sOutput;
            sOutput.enabled = 0; %Computation time-costly, no output, only redraw at the end!

            settings.setValue('forward', forward);            
            [tspan, optionsODE] = obj.buildOptions(settings);
            system = settings.system;
            handles = system.handle();
            x0 = settings.coord;
            x0 = x0(:); %make column vector
            ortho_steps = settings.con_ortho; 
            storage_steps = settings.con_storage;
            param = num2cell(settings.parameters);
            method = str2func(obj.methodname);
            ParamVals=settings.ParamVals;
            % the method variable is actually just the integrator, the LYAP
            % function should be imported here and then method should be
            % passed into that

            % We choose between two modes depending on the parameter
            % settings (A) If there is no active parameter, then we compute
            % Lyapunov Exponents for a single parameter value, but we
            % include the intermediate estimates as the Event so one can
            % check the accuracy of the computation. (B) If there is a
            % single active parameter, we compute Lyapunov Exponents for
            % the entire parameter range and store only the final state and
            % exponents.
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            sOutput.enabled = 0; %The comput            
            if isempty(activeParams)
%              if isempty( optionsODE.Events )
                %In fact, we treat computing Lyapunov Exponents as an Event
                %function as we set a time for re-orthogonalization.
                [t,y,tLE,yLE,tf,yf] = LE(handles, tspan, x0, optionsODE, param, method, ortho_steps, storage_steps);
                % [t, y, Lexp] = LE(length(x0), rhs_ext_maker(x0,handles,param), method, tspan(1), 0.5, tspan(end), x0, 0);
                tE = tf; yE = yf; iE = 1;
                yLE=yLE'; %For correct plotting selection
                % for extend, the LE function should continue from the last
                % lyapunov vectors, not start from unit vectors (as it
                % currently does) 
%              else
%To be adapted, the Lyapunov exponent event is set first, and user functions come next, so this iE-index set one higher.                
%                global Matcont_call
%                Matcont_call = @(tspan, x0, opt) method(handles{2}, tspan, x0, opt, param{:});  
%                [t, y, tE, yE, iE] = method(handles{2}, tspan, x0, optionsODE, param{:});
%              end
              solution = SimCompSolution(settings, obj, tLE, yLE, [], [], [],  method, tspan, x0', optionsODE, cell2mat(param));
              % obj = CompSolution(settings, compbranch);
              % obj.method = method; obj.tspan = tspan; obj.x0 = x0; obj.options = options; obj.param = param;
              %obj.t=t;obj.y=y;
              solution.t=tLE;solution.y=yLE;
              % solution = SimCompSolution(settings, obj, t, y, tLE, yLE, iE,  method, tspan, x0', optionsODE, cell2mat(param));
            elseif length(activeParams)==1
              pdim=length(ParamVals);
              n=length(x0);
              FinalRes=nan(1+2*n,pdim);
              for id=1:pdim
                param{activeParams}=ParamVals(id); 
                [t,y,tLE,yLE,tf,yf] = LE(handles, tspan, x0, optionsODE, param, method, ortho_steps, storage_steps);
                if settings.ReUsePrevComp
                   x0=y(end,:);x0=x0(:);
                end
                FinalRes(:,id)=[ParamVals(id) y(end,:) yLE(end,:)];                
              end
              iE=[];
              solution = SimCompSolution(settings, obj, tLE, yLE, [], [], [],  method, tspan, x0', optionsODE, cell2mat(param));
              solution.t=ParamVals;solution.y=FinalRes(n+2:end,:);
              % disp(FinalRes)
            elseif length(activeParams)>1
              sprintf('The input for at most %s parameter can be an array.', num2str(obj.nrActiveMax))
            end

        end

        function msg = check(obj, settings)
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            if length(activeParams) > obj.nrActiveMax
                msg = sprintf('The input for at most %s parameter can be an array.', num2str(obj.nrActiveMax));
            else
               msg = ''; 
            end
        end



    end 

    methods(Access=protected)
        
        function [tspan, options] = buildOptions(~, settings)
            system = settings.system;
            handles = system.handle();
            
            options = odeset(); %empty
            
            
            options = odeset(options, 'MaxStep', settings.MaxStepSize_sim);   %value is [] when not set in GUI
            options = odeset(options, 'InitialStep', settings.InitStepSize_sim);   %value is [] when not set in GUI
            options = odeset(options, 'RelTol', settings.RelTolerance);
            options = odeset(options, 'AbsTol', settings.AbsTolerance);
            options = odeset(options, 'Refine', settings.Refine);
            options = odeset(options, 'NormControl', CLbool2text(settings.Normcontrol));
            
            options = odeset(options, 'Jacobian', handles{3});
            options = odeset(options, 'Hessians', handles{5});
            
            ef = settings.getSetting('eventfunction');
            if ~isempty(ef) && ef.isVisible() && ~isempty(ef.getValue())
                eventfunction = ef.getValue();
                options = odeset(options, 'Events', eventfunction);
            end
            
            t = settings.time;
            interval = settings.Interval;
            tspan = [t, t+interval];
            if ~settings.forward
                interval = abs(diff(tspan));
                tspan = [tspan(1), tspan(1) - interval];
            end
            
            global sOutput
            if ~isempty(sOutput) && sOutput.enabled
                options = odeset(options, 'OutputFcn', @GUIODEOutput);
                if ~isempty(ef)
                    global Matcont_options
                    Matcont_options = options;
                end
            end
 
        end

        function install(~, settings, name, initval, restrict, catdata)
            settings.addSetting(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingsHelp.getHelp(name) ));
        end

        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            valid = true;
            msg = '';
            if length(activeParams) ~= obj.nrActiveMax
                valid = false;
                msg = sprintf('At most %s parameter can be an array.', num2str(obj.nrActiveMax));
            end
        end


    end


end 
