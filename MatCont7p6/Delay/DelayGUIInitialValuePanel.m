classdef DelayGUIInitialValuePanel

    properties
        panelhandle;
        slider;
        system;
        settings;
    end

    methods
        function obj = DelayGUIInitialValuePanel(session, varargin)
            windowlabel = 'delayinit';
            wm = session.windowmanager;
            parent = wm.demandWindow(windowlabel);
            obj.panelhandle = uipanel(parent);

            obj.system = session.getSystem();
            obj.settings = session.settings;
            delay = obj.system.getDelayInfo();

            dim = size(delay.coordinates,1);
            grid = cell(dim+4, 1);
            options = {'Units', 'Pixels', 'Visible', 'on'};
            suggestions = struct('labelsize', 150);

            grid{1,1} = GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10);
            sectionsettings = DefaultValues.SECTIONNAMESETTING;
            grid{2,1} = uicontrol(obj.panelhandle, 'Style', 'text', 'String', ...
                {'For each coordinate, enter a function for the'; ...
                'initial point. The function can depend on a single'; ...
                'variable, which is evaluated in the delay interval.'; ...
                'All expressions should evaluate to a numerical'; ...
                'value. If the component is constant, the input can'; ...
                'simply be that scalar value. Expressions can be'; ...
                'parameter-dependent. If the delays depend on the'; ...
                'parameters, set them first. To use the values of'; ...
                'the current initial point for a component, leave'; ...
                'the corresponding field empty.'}, ...
                'Max', 3, sectionsettings{:});

            for i = 1:dim
                fieldname = ['co_', delay.coordinates{i,1}, '_init_func'];
                obj.settings.addSetting(fieldname, CLSetting(delay.coordinates{i,1}, '', InputRestrictions.NONE, 2, 20, 1, '~~~~'));
                obj.settings.setVisible(fieldname, false);
                field = obj.settings.getSetting(fieldname);
                grid{i+2,1} = field.renderGUI(session, obj.settings, fieldname, obj.panelhandle, options, suggestions);
                grid{i+2,1}.grid{2}.handle.Callback = @obj.callbackField;
                grid{i+2,1}.grid{2}.handle.Tag = num2str(i);
            end

            grid{end-1,1} = GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10);
            grid{end,1} = uicontrol(obj.panelhandle,'Style','pushbutton','String','OK','Callback',@obj.okCallback);

            box = LayoutBox(grid);
            obj.panelhandle.UserData = box;
            obj.panelhandle.ResizeFcn = @(o, e) obj.onResize(o);

            set(parent, 'Visible' , 'on');
        end

        function obj = onResize(obj, source)
            box = source.UserData;
            delete(obj.slider);
            source.Units = 'pixels';
            obj.slider = box.doSliderLayout(source);
            source.Units = 'normalized';
        end

        function callbackField(obj,src,~)
            % Check that input is a numeric expression or a function handle
            % depending on one variable evaluating to a numeric value.
            % In both cases the expression may depend on the parameters.

            string = src.String;
            index = str2double(src.Tag);
            delay = obj.system.getDelayInfo();
            fieldname = ['co_', delay.coordinates{index,1}, '_init_func'];
            string = renameforsym(string,[],strjoin(obj.system.parameters,','));

            pa = obj.system.parameters;
            for i = 1:length(pa)
                eval( ['par_',pa{i},'=obj.settings.getValue(''pa_',pa{i},''');']);
            end

            field = obj.settings.getSetting(fieldname);
            valid = false;
            errormsg = 'Input must be a numeric expression or a function handle of one variable evaluating to a numeric value, possibly depending on the parameters.';
            try
               string=strtrim(string);
                if isempty(strtrim(string))
                    field.setValue('');
                    obj.settings.refresh();
                else
                    evstring = eval(string);
                    if isnumeric(evstring)
                        string = ['@(t)',string];
                        valid = field.setValue(src.String);
                    elseif isa(evstring,'function_handle')
                        testnum = evstring(0);
                        if isnumeric(testnum)
                            valid = field.setValue(src.String);
                        end
                    end

                    if (~valid)
                        fprintf(2, sprintf('[%s] ERROR(%s): %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), field.displayname, errormsg, string));
                        performErrorDisplay(src);
                        set(src, 'String', field.toString()); %restore original value
                    else
                        obj.settings.refresh();
                    end
                end
            catch error
                fprintf(2, sprintf('[%s] ERROR(%s): %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), field.displayname, error.message, string));
                performErrorDisplay(src);
                set(src, 'String', field.toString()); %restore original value
            end
        end

        function okCallback(obj,~,~)
            % Create new vector of coordinates, update settings and close.

            pa = obj.system.parameters;
            for i = 1:length(pa)
                eval( ['par_',pa{i},'=obj.settings.getValue(''pa_',pa{i},''');']);
            end

            delay = obj.system.getDelayInfo();
            dim = size(delay.coordinates,1);
            fhandles = cell(dim,1);
            for i = 1:dim
                fieldname = ['co_', delay.coordinates{i,1}, '_init_func'];
                field = obj.settings.getSetting(fieldname);
                string = renameforsym(field.value,[],strjoin(obj.system.parameters,','));
                if isempty(string)
                    fhandles{i} = '';
                else
                    fhandles{i} = eval(string);
                    if isnumeric(fhandles{i})
                        fhandles{i} = eval(['@(t)',string]);
                    end
                end
            end
            [~, delay_handles] = feval(obj.system.handle);
            [kmrgd, tau_max] = delay_handles{1}(fhandles, obj.settings.parameters, obj.settings.coord);
            if isnan(tau_max) || tau_max == 0
                error_message = 'With the current parameters, the delay is 0 or there are negative delays. Please set the parameters before you construct the initial point from a function.';
                fprintf(2, sprintf('[%s] ERROR: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), error_message));
                errordlg(error_message, 'Error');
                return
            end
            obj.settings.setValue('coord',kmrgd);

            delete(obj.panelhandle.Parent);
        end
    end
end

function performErrorDisplay(obj)
set(obj, 'BackgroundColor', [1 0.3 0.3]);
pause(0.3)
set(obj, 'BackgroundColor', [1 1 1]);
set(obj, 'BackgroundColor', get(obj, 'BackgroundColor'));
end
