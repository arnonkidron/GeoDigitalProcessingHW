urls = askUserForMultipleOFFfiles();

numMeshHW2es = size(urls,2);

kernelDimFtoV = zeros(numMeshHW2es,1) - 1;
kernelDimVtoF = zeros(numMeshHW2es,1) - 1;
name = cell(numMeshHW2es,1);
V  = zeros(numMeshHW2es,1) - 1;

i = 1;
for url_1 = urls
    url = url_1{1};
    mesh = MeshBasic(url);
    
    [~, name{i}, ~] = fileparts(mesh.Name);
    V(i) = mesh.numV;
    if(mesh.numV >= 10000)
        i = i + 1;
        continue;
    end
    
    mesh = MeshHW2(url);
    
    nullSpace = null(full(mesh.InterpolantFtoV));
    kernelDimFtoV(i) = size(nullSpace,2);
    nullSpace = null(full(mesh.InterpolantVtoF));
    kernelDimVtoF(i) = size(nullSpace,2);
    i = i + 1;
end

table(name, V, kernelDimFtoV, kernelDimVtoF)