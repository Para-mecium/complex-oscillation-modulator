function save_figure(figHandle, fileBase)
figureDir = fileparts(fileBase);
if ~isfolder(figureDir)
    mkdir(figureDir);
end

axesHandles = findall(figHandle, 'Type', 'axes');
for i = 1:numel(axesHandles)
    if isprop(axesHandles(i), 'Toolbar') && ~isempty(axesHandles(i).Toolbar)
        axesHandles(i).Toolbar.Visible = 'off';
    end
end

exportgraphics(figHandle, [fileBase '.png'], 'Resolution', 300);
exportgraphics(figHandle, [fileBase '.pdf'], 'ContentType', 'vector');
end
