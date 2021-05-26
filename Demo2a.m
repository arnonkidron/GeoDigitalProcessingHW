% Gradient
%% valences
for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});

    vertexFunc = GetValences(mesh);
    [~, ~, arrows] = RenderGradient(mesh, vertexFunc);
    title("Gradient of Valences");

    colormap flag
    set(arrows, 'color', 'Yellow');
end

%% distance from center
for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});

    vertexFunc = vecnorm(mesh.Vertices - mean(mesh.Vertices, 1), 2,2);
    [~, ~, arrows] = RenderGradient(mesh, vertexFunc);
    title("Gradient of Distance from Origin");

    colormap flag
    set(arrows, 'color', 'Yellow');
end


    
% Gradient
%% start

for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    
    RenderDivergence(mesh, getVectorFieldUp(mesh));
    title("Vector Field where all vectors point Up");
end


%% normals
for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    
    RenderDivergence(mesh, mesh.FaceNormals);
    title("Vector Field of Face Normals");
    caxis([-1 +1]);
end

%% center to faces
for url = askUserForMultipleOFFfiles()
    mesh = Mesh(url{1});
    
    RenderDivergence(mesh, getCenterToFaces(mesh));
    title("Vector Field of Mesh Center to Faces");
end


%% custom vector field

dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
meshName = "tri_quad_grid.off";


mesh = Mesh(dir + meshName);

vectorField = getCustomVectorField(mesh, ...
    @sin, ...
    @cos, ...
    @tanh ...
    );
RenderDivergence(mesh, vectorField);

%% vector fields


function vectorField = getVectorFieldUp(mesh)
    vectorField = zeros(mesh.numF, 3);
    vectorField(:,3) = 1;
end

function vectorField = getCenterToFaces(mesh)
    faceCenters = getTriangleCenters(mesh);
    center = mean(faceCenters, 1);
    vectorField = faceCenters - center;
end

function vectorField = getCustomVectorField(mesh, f1, f2, f3)
    faceCenters = getTriangleCenters(mesh);
    vectorField = [ ...
        f1(faceCenters(:,1)), ...
        f2(faceCenters(:,2)), ...
        f3(faceCenters(:,3)) ...
    ];
end
