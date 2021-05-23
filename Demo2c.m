mesh = Mesh("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\tri_quad_grid.off");

myRender(mesh, mesh.VertexAreas);
title("Vertex Areas");

[fig, hPatch] = Render(mesh, "White");
view(2);
title("Click on a vertex to compute its row of the Laplacian matrix");

cb = colorbar;
vals = unique(mesh.Laplacian);
vals = full(vals);
set(cb, 'ticks', vals);
set(cb, 'ylim', [min(vals), max(vals)]);
set(get(fig, 'CurrentAxes'), 'Clim', [min(vals), max(vals)]);
colormap parula

% set callback
setappdata(fig,'mesh',mesh);
setappdata(fig,'cb',cb);
set(gca,'ButtonDownFcn',@mouseCB); 
set(hPatch,'ButtonDownFcn',@mouseCB);


function myRender(mesh, f)
    [fig,~] = Render(mesh, f);
    cb = colorbar;
    set(cb, 'ticks', full(unique(f)));
    colormap parula
    view(2);
end

function [] = mouseCB(hObject,event)
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
    
    % replot
    row = mesh.Laplacian(v,:);
    row = full(row)';
    p.FaceVertexCData = row;
    p.FaceColor = 'interp';
    
    % show only the ticks for the current colors
    cb = getappdata(fig, 'cb');
    set(cb, 'ticks', unique(row));
end
