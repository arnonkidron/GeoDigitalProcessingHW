% sometimes we need to compute more than one eigenvalues, so that
% matlab won't complain that the convergence rate is too slow
NUM_EIGENVALUES = 3;

urls = askUserForMultipleOFFfiles();
urls_size = size(urls,2);

names = cell(urls_size, 1);
smallestEigs = zeros(urls_size, NUM_EIGENVALUES);
i = 1;

for url = urls
    mesh = Mesh(url{1});
    W = mesh.CotangentWeights;
    smallestEigs(i,:) = eigs(W, NUM_EIGENVALUES, 'smallestreal');
    
    names{i} = mesh.Name;
    i = i + 1;
end


smallestEigenvalue = smallestEigs(:,1);
table(names, smallestEigenvalue)
