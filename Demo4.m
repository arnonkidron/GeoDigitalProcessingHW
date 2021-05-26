for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    Render(mesh, full(mesh.MeanCurvature));
    colormap hsv
    Render(mesh, full(mesh.GaussianCurvature));
end
