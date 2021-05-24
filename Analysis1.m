urls = askUserForMultipleOFFfiles();
urls_size = size(urls,2);

names = cell(urls_size, 1);

WfNorm = zeros(urls_size, 1);
SymDiffNorm = zeros(urls_size, 1);

numVertices = zeros(urls_size, 1);
numEdges = zeros(urls_size, 1);
numNonZero = zeros(urls_size, 1);
i = 1;

for url = urls
    mesh = Mesh(url{1});
    W = mesh.CotangentWeights;
    W = round(W, 9);
    
    
    % (NULL): show that Wf = 0 for any constant function f
    f = ones(mesh.numV,1);
    WfNorm(i) = norm(sparse(W * f), 'fro');
    
    % (SYM): show that W - W' = 0 
    SymDiffNorm(i) = norm(W - W', 'fro');
    
    % (LOC)
    numNonZero(i) = nnz(W);
    numEdges(i) = mesh.numE;
    numVertices(i) = mesh.numV;
    
    names{i} = mesh.Name;
    i = i + 1;
end

% (NULL) & (SYM)
NULL_SYM_table = table(names, WfNorm, SymDiffNorm)

% (LOC)
numElements = numVertices + 2 * numEdges;
rate = numNonZero ./ numElements;
LOC_table = table(names, numNonZero, numElements, rate)

