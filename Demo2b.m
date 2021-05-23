disp("Frobenius distance between the two laplacian matrices we computed:");


urls = askUserForMultipleOFFfiles();
urls_size = size(urls,2);

names = cell(urls_size, 1);
distance = zeros(urls_size, 1);
i = 1;

for url = urls
    mesh = Mesh(url{1});
    
    cotLaplacian = LaplacianByCotFormula(mesh);
    distance(i) = norm(mesh.Laplacian - cotLaplacian, 'fro');
    
    names{i} = mesh.Name;
    i = i + 1;
end

table(names, distance)