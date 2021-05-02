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
    analyze4(filename{1});
end


function fig = analyze4(filename)
    mesh = MeshWithoutAreaStuff(filename);

    A = mesh.Adjacency - mesh.Adjacency';
    isVBoundary = any(A ~= 0, 2);
    isVnotBoundary = not(isVBoundary);
    
    if(all(isVnotBoundary))
        [~,name,~] = fileparts(filename);
        disp(name + " has no boundary");
        return;
    end
    
    isVnotBoundary = full(isVnotBoundary * 1);

    fig = RenderWireframe(mesh, isVnotBoundary);
    set(colorbar, 'ticks', min(isVnotBoundary):max(isVnotBoundary));
    colormap bone;
    
     [~,name,~] = fileparts(filename);
     set(fig, 'Name', name);
end

