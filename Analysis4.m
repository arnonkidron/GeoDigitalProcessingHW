urls = askUserForMultipleOFFfiles();

for url = urls
    analyze4(url{1});
end

function fig = analyze4(url)
    mesh = MeshWithoutAreaStuff(url);

    A = mesh.Adjacency - mesh.Adjacency';
    isVBoundary = any(A ~= 0, 2);
    isVnotBoundary = not(isVBoundary);
    
    if(all(isVnotBoundary))
        [~,name,~] = fileparts(url);
        disp(name + " has no boundary");
        return;
    end
    
    isVnotBoundary = full(isVnotBoundary * 1);

    fig = RenderWireframe(mesh, isVnotBoundary);
    set(colorbar, 'ticks', min(isVnotBoundary):max(isVnotBoundary));
    colormap bone;
    
     [~,name,~] = fileparts(url);
     set(fig, 'Name', name);
end

