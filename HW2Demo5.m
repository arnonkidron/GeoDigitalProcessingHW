urls = askUserForMultipleOFFfiles();

for url = urls
    [~,fileName,format] = fileparts(url{1});
    mesh = MeshHW2(fileName + ".off");
    extrema = [min([mesh.TriangleAreas; mesh.VertexAreas]), max([mesh.TriangleAreas; mesh.VertexAreas])];
    myRender(mesh, mesh.TriangleAreas, "Triangle Areas", extrema);
    myRender(mesh, mesh.VertexAreas, "Vertex Areas", extrema);
end


function myRender(mesh, f, titleText, y_limits)
    fig = Render(mesh, f);
    set(get(fig, 'CurrentAxes'), 'Clim', y_limits);
    cb = colorbar;
    set(cb, 'ticks', unique(f));
    set(cb, 'ylim', y_limits);
    colormap parula
    
    title(titleText);
end