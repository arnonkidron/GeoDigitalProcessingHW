urls = askUserForMultipleOFFfiles();

for url = urls
    analyze2(url{1});
end

function fig = analyze2(url)
    mesh = MeshWithoutAreaStuff(url);

    A = mesh.Adjacency + mesh.Adjacency';
    valence = full(sum(A ~= 0, 2));

    fig = Render(mesh, valence);
    set(colorbar, 'ticks', min(valence):max(valence));
    
    [~,name,~] = fileparts(url);
    set(fig, 'Name', name);
end
