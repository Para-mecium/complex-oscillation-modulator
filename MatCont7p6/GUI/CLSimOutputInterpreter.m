%default output interpreter for continuation curves
classdef CLSimOutputInterpreter
    
    methods
        
        %{
        %unused.
        function [coord, param, period, meta] = getPoint(~, solution, index)
        end
        %}
        
        function plotter = getPlotter(~)
           plotter = @GUIPlotOutputterSim;
        end
        
        
        function renderer = getPreviewPanelRenderer(~)
           renderer = @(sol, pp, pm) GUIPreviewRendererSolutions.renderSimulation(sol, pp, pm);
        end
        
        function [omap, numconfig, plotsel] = interpret(~, solution, settings, compbranch)
            omap = []; %used for continuation curves.

            numconfig = settings.numericconfig;
            numconfig.reset();
            plotsel = PlotSelections();
            
            system = settings.system;
            coordinates = system.getCoordinates();
            dim = length(coordinates);
            parameters = system.getParameters();
            parametervalues = settings.getValue('parameters');
            
          if ~isa(compbranch,'SimConf_LyapExp')
            numconfig.declareCategory('coordinates', 'Coordinates',1e0, true);
            plotsel.declareCategory('coordinates', 'Coordinates', 1e0);
            
            labels = cell(dim, 2);
            for index = 1:dim
                labels(index, :) = {coordinates{index}, @(t, y, s, ind, i) y(i, index)};
                plotsel.declareSubCategory('coordinates', coordinates{index}, @(t, y, s, ind, i) y(ind, index))
            end
            numconfig.setLabels('coordinates', labels);            

           % assign system to variable, otherwise methods cannot be
            % used (see CLSettings>subsref, line 159: why is going to
            % levels other than 1 and 3 forbidden?)
            system = settings.getVal('system', []);
            if ~isempty(system)
                delay = system.getDelayInfo;
                if isempty(delay)
                    % not an imported delay equation
                    hasrenewal = 0;
                else
                    hasrenewal = any(cell2mat(delay.coordinates(:,3))==0);
                end
                if hasrenewal
                    [~, delay_handles] = feval(system.handle);
                    % for backward compatibility with 7p5 uncomment all
                    % of the following if block
                    % if isa(delay_handles, 'cell')
                        delay.rhs_r = delay_handles{2};
                    % elseif isa(delay_handles, 'function_handle')
                    %     delay.rhs_r = delay_handles;
                    % end
                    % indices of coordinates described by renewal equations
                    delay.ind_r = find(cell2mat(delay.coordinates(:,3))==0);
                    % dimension of the renewal part of the system
                    delay.dim_r = length(delay.ind_r);

                    numconfig.declareCategory('reconstructed', 'Reconstructed coordinates (renewal equations)',1e0, true);
                    plotsel.declareCategory('reconstructed', 'Reconstructed coordinates (renewal equations)', 1e0);

                    labels = cell(delay.dim_r, 2);
                    for index = 1:delay.dim_r
                        labels(index, :) = {delay.coordinates{delay.ind_r(index),1}, @(t, y, s, ind, i) (delay.rhs_r((y(i, :))',parametervalues,[],delay.ind_r(index),0,0))'};
                        plotsel.declareSubCategory('reconstructed', delay.coordinates{delay.ind_r(index),1}, @(t, y, s, ind, i) (delay.rhs_r((y(ind, :))',parametervalues,[],delay.ind_r(index),0,0))')
                    end
                    numconfig.setLabels('reconstructed', labels);
                end
            end
        end


            numconfig.declareCategory('parameters', 'Parameters',1e2, false)
            plotsel.declareCategory('parameters', 'Parameters',1e2)
            labels = cell(length(parameters), 2);
            
            partemp=settings.getSetting('parameters');
            ap=partemp.getActive();
            for index = 1:length(parameters)
                labels(index, :) = {parameters{index}, @(t, y, s, ind, i) parametervalues(index) };
                if (index ~= ap)
                  plotsel.declareSubCategory('parameters', parameters{index},  @(t, y, s, ind, i) repmat(parametervalues(index), 1, length(ind)));
                else
                  plotsel.declareSubCategory('parameters', parameters{index},  @(t, y, s, ind, i) t(ind));
                end
            end
            numconfig.setLabels('parameters', labels);
            
            timevar = system.getTimeName();
            numconfig.declareCategory('time', 'Time',101, true);
            numconfig.setLabels('time', {timevar,  @(t, y, s, ind, i) t(i)});
            plotsel.declareCategory('time', 'Time', 101);
            plotsel.declareSubCategory('time', timevar, @(t, y, s, ind, i) t(ind));
            
            if isa(compbranch,'SimConf_LyapExp')
                % Numeric Panel is too slow...
                % numconfig.declareCategory('LyapExps','Lyapunov exponents', 2e2);
                plotsel.declareCategory('LyapExp', 'Lyapunov exponents', 2e2); 
                for index = 1:dim
                    labels(index, :) = {strcat('LyapExp',num2str(index)), @(t, y, s, ind, i) y(i, index)};
                    plotsel.declareSubCategory('LyapExp',strcat('LE_',num2str(index')), @(t, y, s, ind, i) y((ind-1)*dim+index));
                end
                % numconfig.setLabels('LyapExps', labels);
                plotsel.declareSubCategory('LyapExp',strcat('LEsAll'), @(t, y, s, ind, i) y);
            end

            if compbranch.hasTarget()
                targetsetting = settings.getSetting('coord_target');
                if ~isempty(targetsetting) && targetsetting.isVisible()
                    target = targetsetting.getValue();
                else
                    target = settings.coord;
                end
                target = target(:)';
                name = compbranch.hasTarget();
                numconfig.declareCategory(name, name,1e4, true);
                numconfig.setLabels(name, {name,  @(t, y, s, ind, i) norm(y(i,:) - target)});
            end
            
            ufdata = system.ufdata;
            if ~isempty(ufdata) && ~isempty(find([ufdata.valid]))
                paramcell = num2cell(parametervalues);
                handles = system.handle(); %evaluate handles
                %out{1} = @init;out{2} = @fun_eval;out{3} = @jacobian;out{4} = @jacobiap;out{5} = @hessians;%
                %out{6} = @hessiansp;out{7} = @der3;out{8} = [];out{9} = [];
                %out{10}= @Userfun1;out{11}= @Userfun2; ... ...             
                handles(1:9) = [];
                numconfig.declareCategory('uf', 'User Functions',1e6, false)
                plotsel.declareCategory('uf', 'User Functions',1e6)
                actives = find([ufdata.valid]);
                labels = cell(length(actives), 2);
                j = 1;
                
                for k = 1:length(ufdata)
                    index = index + 1;
                    if any(k == actives)
                        userfun = handles{k};
                        labels(j, :) = {strip(ufdata(k).label), @(t, y, s, ind, i) userfun(t(i), y(i, :)', paramcell{:})};
                        j = j + 1;
                        plotsel.declareSubCategory('uf',strip(ufdata(k).label) , @(t, y, s, ind, i)  arrayfun(@(z)userfun(t(z), y(z, :)', paramcell{:}), ind));
                    end
                end
                
                numconfig.setLabels('uf', labels);
            end
           
            plotsel.setSolutionLabel(compbranch);
            numconfig.configDone();
        end
    end
end
