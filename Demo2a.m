%%
mesh = Mesh("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\disk.off");
 
vertexFunc = GetValences(mesh);
[~, ~, arrows] = RenderGradient(mesh, vertexFunc);

view(2)
colormap flag
set(arrows, 'color', 'Yellow');

% 
% 
% 
% %mesh = Mesh("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\triangle.off");
% %vertexFunc = [0 0 1]';
% % simple = OldGradient(mesh, vertexFunc);
% % complex = reshape(Gradient(mesh, vertexFunc), [mesh.numF, 3]);
% % RenderVectorField(mesh, vertexFunc, simple);
% % view(2);
% % colormap flag;
% % RenderVectorField(mesh, vertexFunc, complex);
% % view(2);
% % colormap flag;
% 
% 
% lap = Laplacian(mesh, vertexFunc);
% Render(mesh, vertexFunc);
% Render(mesh, lap);

% 
% mesh = Mesh("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\tri_quad_grid.off");
% Render(mesh, []);

%%
% mesh = Mesh([]);
% diff = mesh.laplacian - mesh.cotLaplacian;
% sparse_diff = sparse(diff);

%% start
urls = askUserForMultipleOFFfiles();

for url = urls
    mesh = Mesh(url{1});
    divergenceUp(mesh);
end

function divergenceUp(mesh)
    vectorField = zeros(mesh.numF, 3);
    vectorField(:,3) = 1;
    [fig, p, arrows] = RenderVectorField(mesh, Divergence(mesh, vectorField), vectorField);
    title("Vector Field Up");
end




% More ideas:
% normals
% center to point
% only up
% functions on grid


