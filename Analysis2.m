dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
linkList = [
    "cat.off"
    ];
% TODO: iterate over all files in directory
% ask use for directory

for link = linkList
    mesh = Mesh(dir + link);


    A = mesh.Adjacency + mesh.Adjacency';
    valence = full(sum(A ~= 0, 2));

    Render(mesh, valence);
    set(colorbar, 'ticks', min(valence):max(valence));
end