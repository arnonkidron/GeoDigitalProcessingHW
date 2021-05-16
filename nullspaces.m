urls = askUserForMultipleOFFfiles();

numMeshes = size(urls,2);

kernelDimFtoV = zeros(numMeshes,1) - 1;
kernelDimVtoF = zeros(numMeshes,1) - 1;
name = cell(numMeshes,1);
V  = zeros(numMeshes,1) - 1;

i = 1;
for url_1 = urls
    url = url_1{1};
    mesh = MeshWithoutAreaStuff(url);
    
    [~, name{i}, ~] = fileparts(mesh.Name);
    V(i) = mesh.numV;
    if(mesh.numV >= 10000)
        i = i + 1;
        continue;
    end
    
    mesh = Mesh(url);
    
    nullSpace = null(full(mesh.InterpolantFtoV));
    kernelDimFtoV(i) = size(nullSpace,2);
    nullSpace = null(full(mesh.InterpolantVtoF));
    kernelDimVtoF(i) = size(nullSpace,2);
    i = i + 1;
end

table(name, V, kernelDimFtoV, kernelDimVtoF)