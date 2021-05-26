global EDGE_ALPHA;
EDGE_ALPHA = 0.15;

for url = askUserForMultipleOFFfiles()
    k=10:20:310;
    
    meshSmoother = MeshSmoother(url{1}, max(k));
    
    % render some eigenfunctions
    RenderSomeEigenfunctions(meshSmoother, 9, EDGE_ALPHA);
    
    % do the experiment with the hat function
    f = GetHatFunction(meshSmoother, 255);
    experiment(meshSmoother, f, k);
    
    % do the experiment with the maximal eigenfunction
    f = GetMaximalEigenfunction(meshSmoother);
    experiment(meshSmoother, f, k);
end

function experiment(meshSmoother, f, k)
    global EDGE_ALPHA;
    size_k = size(k,2);
    
    errors = zeros(size_k,1);
    gs = zeros(meshSmoother.Mesh.numV,size_k);
    titles = cell(size_k,1);
    
    for i=1:size_k
        ki = k(i);
        [g, error] = ProjectAndBack(meshSmoother, f, ki);
        errors(i) =  error;
        gs(:,i) = g;
        titles{i} = sprintf("k = %d", ki);
    end
    
    indices_to_render = [1:3:16];
    RenderSeveralFunctions(meshSmoother.Mesh, ...
        [gs(:, indices_to_render) f], ...
        {titles{indices_to_render} "f"}, ...
        EDGE_ALPHA);
    
    figure
    plot(k, errors, 'o-');
    xlabel("k");
    ylabel("error");
    
end


