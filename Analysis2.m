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
    analyze2(filename{1});
end

function analyze2(filename)
    mesh = MeshWithoutAreaStuff(filename);

    A = mesh.Adjacency + mesh.Adjacency';
    valence = full(sum(A ~= 0, 2));

    Render(mesh, valence);
    set(colorbar, 'ticks', min(valence):max(valence));
end
