
for url = askUserForMultipleOFFfiles()
    k=10:20:310;
    
    meshPCA = MeshPCA(url{1}, max(k));
%     Render9Eigenfunctions(meshPCA);

    f = GetHatFunction(meshPCA, 255);
    experiment(meshPCA, f, k);
    
%     f = GetMaximalEigenfunction(meshPCA);
%     experiment(meshPCA, f, k);
end

function experiment(meshPCA, f, k)
    size_k = size(k,2);
    
    errors = zeros(size_k,1);
    gs = zeros(meshPCA.Mesh.numV,size_k);
    
    for i=1:size_k
        ki = k(i);
        [g, error] = ProjectAndBack(meshPCA, f, ki);
        errors(i) =  error;
        gs(:,i) = g;
    end
    
    gs_to_render = [1 4 7 10 13 16];
    RenderSeveralFunctions(meshPCA.Mesh, [gs(:, gs_to_render) f]);
    
    figure
    plot(k, errors, 'o-');
    xlabel("k");
    ylabel("error");
    
end


