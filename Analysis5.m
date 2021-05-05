urls = askUserForMultipleOFFfiles();

for url = urls
    analyze5(url{1});
end

function fig = analyze5(url)
    mesh = MeshWithoutAreaStuff(url);

    [fig, hPatch] = RenderWireframe(mesh, 'Black');
    
    setappdata(fig,'mesh',mesh);
    set(gca,'ButtonDownFcn',@mouseCB5); 
    set(hPatch,'ButtonDownFcn',@mouseCB5);
    
    set(fig, 'Name', get(fig, 'Name') + " - distances from point");
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
