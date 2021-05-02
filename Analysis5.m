% ask the user to select the files to process. He may choose multiple files
[filenames, paths] = uigetfile('*.off', ...
    'Select files to process', ...
    'MultiSelect', 'on');
filenames = fullfile(paths, filenames);

% if the user has selected one file, put it in a cell, so that the code
% afterwards will be the same
if(~iscell(filenames))
    filenames = {filenames};
end

% analyze all of the files
for filename = filenames
    analyze5(filename{1});
end


function fig = analyze5(filename)
    mesh = MeshWithoutAreaStuff(filename);

    [fig, hPatch] = RenderWireframe(mesh, 'Black');
    
    setappdata(fig,'mesh',mesh);
    set(gca,'ButtonDownFcn',@mouseCB5); 
    set(hPatch,'ButtonDownFcn',@mouseCB5);
    
    [~,name,~] = fileparts(filename);
    set(fig, 'Name', name);
end

function [] = mouseCB5(hObject,event)
    % get XYZ coordinates & mesh
    type = get(hObject,'Type');
    if strcmpi(type,'axes')
        axesObject = hObject;
        fprintf('No patch!\n');
    elseif strcmpi(type,'patch')
        axesObject = get(hObject,'Parent');
        p = hObject;
    else
        fprintf('unhandled object type!\n');
        return;
    end
    XYZ = get(axesObject,'CurrentPoint');
    fig = get(axesObject,'Parent');
    mesh = getappdata(fig, 'mesh');
    
%     xPosition = XYZ(1,1);
%     yPosition = XYZ(1,2);
%     zPosition = XYZ(1,3);
%     fprintf('x=%f Y=%f Z=%f \n',xPosition,yPosition,zPosition)

    distanceFromClick = sqrt(sum(bsxfun(@minus, mesh.Vertices, XYZ(1,:)).^2,2));
    [~,v] = min(distanceFromClick);
    
    G = GetWeightedGraph(mesh);
    
    distFromV = distances(G, v)';
    
    % replot
    p.FaceVertexCData = distFromV;
    p.EdgeColor = 'interp';
    
    
    

end
