for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
%     Render(mesh, full(mesh.MeanCurvature));
    colormap hsv
    Render(mesh, full(mesh.GaussianCurvature));
end

%%
mesh = Mesh([]);
%%
v = unique(unique(full(mesh.GaussianCurvature)))

%%
v = sort(v')
%%
u = round(v,4)
u = unique(u)

%%
u'