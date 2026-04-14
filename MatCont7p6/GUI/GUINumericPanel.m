classdef GUINumericPanel < handle
    properties

        panel;
        model;
        layouthandle;
        optionshandle;
        setters = {};
        
        %options =  {'Units', 'Pixels', 'Visible', 'on'};
        %suggestions = struct('labelsize', 150 );
        %widthmargin = 10;
        %session = [];
        settings = [];
        
        eventlistener = [];
        eventlistener2 = [];
        
        slider;
        
    end
      
    methods
        function obj = GUINumericPanel(parent, numericmodel, settings)
            obj.panel = uipanel(parent, 'DeleteFcn' , @(o,e) obj.destructor(), 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            obj.layouthandle = uimenu(parent, 'Label', 'Layout');
            obj.optionshandle = uimenu(parent, 'Label', 'Options');
            obj.settings = settings;
            obj.setupEigMultOptions();
            uimenu(obj.optionshandle, 'Label', 'Eigenvalues/Multipliers', 'Callback', @(src,ev) obj.eigMultOptionsDialog());
            obj.changeNumericModel(numericmodel);
            obj.slider = [];
            
            scrollAmmount = DefaultValues.LETTERDIMENSION(2) * 1;
            set(parent, 'WindowScrollWheelFcn', @(o, e) onScrollEvent(obj, e.VerticalScrollCount, scrollAmmount));

        end
        
        function changeNumericModel(obj, numericmodel)
            obj.model = numericmodel;
            delete(obj.eventlistener);
            delete(obj.eventlistener2);
            obj.eventlistener = obj.model.addlistener('visibilityChanged' , @(srv,ev) obj.setupPanel());
            obj.eventlistener2 = obj.model.addlistener('categoriesChanged' , @(srv,ev) obj.setup());
            obj.setup();
        end
        
        function destructor(obj)
           delete(obj.layouthandle);
           delete(obj.eventlistener);
           delete(obj.eventlistener2);
           delete(obj);
           
        end
        
        
        function setup(obj)
            obj.setupPanel();
            obj.setupLayoutMenu();
        end
        
        function setupLayoutMenu(obj)
            delete(allchild(obj.layouthandle));
            catnames = obj.model.getCategories();
            for i = 1:length(catnames)
                GUISwitchMenuItem(obj.layouthandle, obj.model.getFullName(catnames{i}), @(b) obj.model.setCatVisible(catnames{i}, b), @() obj.model.isCatVisible(catnames{i}), obj.model, 'visibilityChanged');
            end
            
        end
        
        function setupEigMultOptions(obj)
            obj.settings.addSetting('EigMultShow', CLSetting('EigMultShow', DefaultValues.NUMERICOPTIONS.EigMultShow, InputRestrictionCategory({'all','some'}), 0, 6, 1, CLSettingsHelp.getHelp('EigMultShow')));
            obj.settings.addSetting('EigMultShowNum', CLSetting('EigMultShowNum', DefaultValues.NUMERICOPTIONS.EigMultShowNum, InputRestrictions.INT_g0, 0, 6, 2, CLSettingsHelp.getHelp('EigMultShowNum')));
        end

        function eigMultOptionsDialog(obj)
            EigMultShow = obj.settings.getSetting('EigMultShow');
            EigMultShowNum = obj.settings.getSetting('EigMultShowNum');
            d = dialog('Name','Options for Eigenvalues and Multipliers',...
                'Position',[300 300 410 140]);
            bg = uibuttongroup(d);
            uicontrol(bg,'Style','radiobutton',...
                'Position',[20 100 260 20],...
                'String','Show all eigenvalues/multipliers',...
                'Tag','all',...
                'Callback',@callback,...
                'Value',strcmp(EigMultShow.getValue,'all'));
            uicontrol(bg,'Style','radiobutton',...
                'Position',[20 60 90 20],...
                'String','Show only',...
                'Tag','some',...
                'Callback',@callback,...
                'Value',strcmp(EigMultShow.getValue,'some'));
            N = uicontrol(d,'Style','edit',...
                'Position',[120 60 40 20],...
                'Callback',@callback,...
                'String',sprintf('%d',EigMultShowNum.getValue));
            uicontrol(d,'Style','text',...
                'Position',[170 60 220 20],...
                'String','dominant eigenvalues/multipliers');
            uicontrol(d,'Style','pushbutton',...
                'Position',[155 20 100 20],...
                'String','OK',...
                'Callback','delete(gcf)');
            uiwait(d);
            obj.setupPanel();

            function callback(~,~)
                obj.settings.setValue('EigMultShow',bg.SelectedObject.Tag);
                obj.settings.setValue('EigMultShowNum',str2double(N.String));
            end
        end

        function setupPanel(obj)
            delete(allchild(obj.panel));
          
            obj.setters = {};
            catnames = obj.model.getCategories();
            
            grid = cell(0, 2);
            for i = 1:length(catnames)
                if obj.model.isCatVisible(catnames{i})
                    labels = obj.model.getLabels(catnames{i});
                    grid{end+1,2} = GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10);
                    grid{end+1,2} = uicontrol(obj.panel, 'style', 'text', 'String' , obj.model.getFullName(catnames{i}), 'horizontalalignment' , 'center', 'FontAngle', 'italic', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                    
                    if strcmp(catnames{i},'eig') || strcmp(catnames{i},'mult')
                        num_eig = size(labels, 1)/2;
                        EigMultShow = obj.settings.getSetting('EigMultShow');
                        EigMultShowNum = obj.settings.getSetting('EigMultShowNum');
                        if strcmp(EigMultShow.getValue,'some')
                            last_eig = EigMultShowNum.getValue;
                        else
                            last_eig = num_eig;
                        end
                        labels_indices = [1:last_eig,(num_eig+1):(num_eig+last_eig)];
                    else
                        labels_indices = 1:size(labels, 1);
                    end
                    for j = labels_indices
                        datalabel = uicontrol(obj.panel, 'style', 'text', 'String' , '-', 'horizontalalignment' , 'center', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                        obj.setters{end+1} = @(data,s,ind,i) set(datalabel, 'String', sprintf('%.15g', labels{j, 2}(data{:},s,ind,i)));
                        grid(end+1, :) = {uicontrol(obj.panel, 'style', 'text', 'String' , labels{j, 1}, 'horizontalalignment' , 'left', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR) , datalabel};
                    end
                end
                
            end
            mainbox = LayoutBox(grid);
            obj.panel.UserData = mainbox;
            obj.panel.ResizeFcn = @(o, e) obj.onResize(o);
            obj.onResize(obj.panel);
        end
        
        function output(obj, data, s, ind)
            for i = 1:length(obj.setters)
                obj.setters{i}(data, s, ind, ind(end));
            end
        end
        
        function outputPoint(~, varargin)  %not implemented for Numeric.
        end    
        
        function b = isValid(obj)
            b = 1; %check if anything selected?
        end
        function onResize(obj, source)
            box = source.UserData;
            source.Units = 'pixels';
            delete(obj.slider);
            obj.slider = box.doSliderLayout(source);
            source.Units = 'normalized';
        end        
    end

end

function onScrollEvent(panel, direction, amount)
    %direction: 1 if scroll down, -1 if scroll up.
    if isempty(panel.slider); return; end
    
    newvalue = panel.slider.Value - direction*amount;
    newvalue = min(max(panel.slider.Min, newvalue), panel.slider.Max);
    panel.slider.Value = newvalue;
    panel.slider.Callback(panel.slider, []);
 
end
