for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    Render(mesh, full(mesh.MeanCurvature));
%  %   Render(mesh, full(mesh.GaussianCurvature));
end

