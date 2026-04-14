classdef DelayGUICommands

    methods(Static)
        
        function call_new_system(session)
            system = DelayGUICommands.gui_loader(@DelayGUICommands.new);
            if ~isempty(system)
                session.changeSystem(system);
            end
            
        end
        
        function call_edit_system(session)
            
            system = DelayGUICommands.gui_loader(@() DelayGUICommands.edit(session.getSystem().getName()));
            if ~isempty(system)
                session.changeSystem(system);
            end
            
        end
        
        function sys = gui_loader(systemcmd)
            global gds;
            gds = [];
            fhandle = systemcmd();
       
            if ~isempty(fhandle) && isvalid(fhandle); uiwait(fhandle); end
            
            if ~isempty(gds) && isfield(gds, 'ok') && gds.ok
                sysname = gds.system;
                [path_sys, ~, ~] = fileparts(which('standard.m'));
                syspath = fullfile(path_sys, [sysname '.mat']);
                sys = CLSystem(syspath);
                
            else
                sys = [];
            end
            
        end
        
        function outhandle = new()
            global path_sys
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            import_delay('init');
            fhandle = import_delay();
            if nargout == 1
               outhandle = fhandle; 
            end
        end
        
        function outhandle = edit(systemname)
            
            if (~ischar(systemname))
                systemname = func2str(systemname);
            end
            
            global path_sys;
            global gds;
            [path_sys, ~, ~] = fileparts(which('standard.m'));
            path_sys = [path_sys filesep()];
            load( [path_sys  systemname '.mat' ]);  %overwrites gds.
            
            if ~isempty(gds) && ~isfield(gds, 'delay')
                errordlg('The system is not an imported delay equation: you should edit it with MatCont''s standard edit window.','Error')
                outhandle = [];
                gds = [];
                return
            end
            
            fhandle = import_delay();
            if nargout == 1
                outhandle = fhandle;
            end
        end
        
    end
end
