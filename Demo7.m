urls = askUserForMultipleOFFfiles();
urls_size = size(urls,2);

names = cell(urls_size, 1);
boundary_edges = zeros(urls_size, 1);
genus = zeros(urls_size, 1);
x = zeros(urls_size, 1);

i = 1;

for url = urls
    [~,fileName,format] = fileparts(url{1});
    mesh = MeshWithoutAreaStuff(fileName + ".off");
    
    names{i} = fileName;
    boundary_edges(i) = CalcBoundaryEdges(mesh);
    genus(i) = CalcGenus(mesh);
    x(i) =  mesh.numV - mesh.numE + mesh.numF;
    i = i+1;
end

table(names,boundary_edges,genus,x)

