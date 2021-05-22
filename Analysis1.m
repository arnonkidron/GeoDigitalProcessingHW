urls = askUserForMultipleOFFfiles();

for url = urls
    mesh = Mesh(url{1});
    W = mesh.CotangentWeights;
    % (NULL): show that Wf = 0 for any constant function f
    f = ones(mesh.numV,1);
    disp("(NULL): the result should be a zero vector");
    disp(sparse(W * f));
    
    disp("(SYM): the result should be zero");
    disp(norm(W - W', 'fro'));
end
