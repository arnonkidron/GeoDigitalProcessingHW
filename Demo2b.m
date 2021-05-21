mesh = Mesh([]);
diff = mesh.laplacian - mesh.cotLaplacian;
sparse_diff = sparse(diff);

isAntiSymmetric(diff)

function check = isSymmetric(A)
    check = A == A';
end
