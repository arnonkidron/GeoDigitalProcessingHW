for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    
    [fig, hPatch] = Render(mesh, "White");
    view(2);
    title("Click on a vertex to compute its row of the W matrix");

    cb = colorbar;
    vals = unique(mesh.CotangentWeights);
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
    row = mesh.CotangentWeights(v,:);
    row = full(row)';
    p.FaceVertexCData = row;
    p.FaceColor = 'interp';
    
    % show only the ticks for the current colors
    cb = getappdata(fig, 'cb');
    set(cb, 'ticks', unique(row));
end




%% print the names of meshes that have positive edge weights, or negative vertex weights, e.g. the parallelogram with obtuse angles
for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    W = mesh.CotangentWeights;
    withoutDiagonal = full(W) - diag(diag(W));
    negativeOnDiagonal = nnz(diag(W) < 0 - 1e-9);
    positiveOffDiagonal = nnz(withoutDiagonal > 0 + 1e-9);


    if (negativeOnDiagonal + positiveOffDiagonal > 0)
        disp(mesh.Name);
    end
end



