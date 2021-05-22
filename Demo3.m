for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    normals = mesh.VertexNormals;
    RenderVectorField(mesh, [], normals);
end


