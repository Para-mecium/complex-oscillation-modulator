classdef CLSettingCoordinates < CLSetting

    properties
       coordinates 
    end
    
    methods
        function obj = CLSettingCoordinates(settings, coordinates)
            dim = length(coordinates);
            obj = obj@CLSetting('coord', zeros(1, dim), InputRestrictions.vector(dim), 2, 1, length(coordinates) + 1, '~~~');
            obj.coordinates = coordinates;
            obj.displayname = sprintf('alternative\ncoordinate array');

            if ~isempty(settings)
                obj.installGhosts(settings);
            end
        end
        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            grid = cell(2, 1);

            % initialize EditBox while the value of 'coord' is a vector of
            % zeros in order to avoid an overly wide box
            coord_value = settings.getValue('coord');
            settings.setValue('coord', zeros(size(coord_value)));
            grid{1} = renderGUI@CLSetting(obj, session, settings, label, panelhandle, options, suggestions);
            settings.setValue('coord',coord_value);

            setup_button_visible = false;
            setup_button_enable = 'off';
            % assign system to variable, otherwise methods cannot
            % be used (see CLSettings>subsref, line 159: why is going
            % to levels other than 1 and 3 forbidden?)
            system = settings.getVal('system', []);
            if ~isempty(system)
                delay = system.getDelayInfo;
                if ~isempty(delay)
                    % imported delay equation
                    setup_button_visible = true;
                    hasrenewal = any(cell2mat(delay.coordinates(:,3))==0);
                    if nargout(system.handle)>1
                        [~, delay_handles] = feval(system.handle);
                        if isa(delay_handles, 'cell') && (~isscalar(delay_handles) || ~hasrenewal)
                            setup_button_enable = 'on';
                        end
                    end
                end
            end
            grid{2} = uicontrol(panelhandle, 'Style', 'pushbutton', ...
                'Visible', setup_button_visible, 'Enable', setup_button_enable, ...
                'String', 'Setup initial point from function', ...
                'Callback', @(o,e) DelayGUIInitialValuePanel(session));
            box = LayoutBox(grid);  
        end
        
        function installGhosts(obj, settings)
            for i = 1:length(obj.coordinates)
                coordinate = obj.coordinates{i};
                settings.addSetting(['co_', coordinate], CLSettingCoordinate(obj, i));
            end
        end
        
       function newobj = copy(obj, newsettings)
           newobj = CLSettingCoordinates(newsettings, obj.coordinates);
           newobj.value = obj.value;
           obj.copyOver(newobj);
       end        
        
       function b = sanityCheck(obj, settings)
            b = obj.internalSanityCheck(settings);
           
            if ~b
                for i = 1:length(obj.coordinates)
                    coordinate = obj.coordinates{i};
                    settings.removeSetting(['co_', coordinate]);
                end
            end
            
       end
       
       function b = internalSanityCheck(obj, settings)
          system = settings.system;
          if isempty(system); b = 0; return; end
          coordinates = system.getCoordinates();
          if length(coordinates) ~= length(obj.coordinates); b = 0; return; end
          
          for k = 1:length(obj.coordinates)
                if ~strcmp(obj.coordinates{k}, coordinates{k})
                    b = 0;
                    return
                end
          end
          b = 1;
       end          
    end
end
