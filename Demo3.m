urls = askUserForMultipleOFFfiles();

for url = urls
    mesh = Mesh(url{1});
    normals = ComputeVertexNormals(mesh);
    RenderVectorField(mesh, [], normals);
end

