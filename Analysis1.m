% collect data
dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
avgEdgeLengths = zeros(6,1);
avgCenterErrors = zeros(6,1);
for i = 0:5
  %  fprintf("sphere %d\n", i);
    mesh = MeshWithoutAreaStuff(dir + "sphere_s" + i + ".off");
    avgEdgeLengths(i + 1) = CalcAvgEdgeLength(mesh);
    avgCenterErrors(i + 1) = CalcAvgCenterError(mesh);
    
    [~, p] = Render(mesh, vecnorm(mesh.Vertices,2,2));
    p.EdgeAlpha = 0.2 - i / 30;
end

convergenceQuotient = zeros(5,1);
for i = 1:5
    convergenceQuotient(i) = avgCenterErrors(i+1) / avgCenterErrors(i);
end
disp(convergenceQuotient);

% plot
figure
p = plot(avgEdgeLengths, avgCenterErrors, 'o- red', 'LineWidth', 2);
xlabel("Average Edge Length of Mesh");
ylabel("Average Triangulation Error of Mesh");
title("Triangulation Error");


% functions

function center = getTriangleCenter(mesh, f)
    % input: a mesh
    %        a face index
    % output: its center point
    index1 = mesh.Faces(f, 1);
    index2 = mesh.Faces(f, 2);
    index3 = mesh.Faces(f, 3);
    
    v1 = mesh.Vertices(index1, :);
    v2 = mesh.Vertices(index2, :);
    v3 = mesh.Vertices(index3, :);
    
    center = mean([v1; v2; v3], 1);
end

function projection = projectToUnitSphere(p)
    % input: a point
    % output: its projection onto the unit sphere
    projection = p / norm(p);
end


function avgCenterError = CalcAvgCenterError(mesh)
    % input: a mesh that is supposed to approximate the unit sphere
    % output: the average distance between its face centers and their
    % projections onto the unit sphere

    totalError = 0;
    for f = 1:mesh.numF
        center = getTriangleCenter(mesh, f);
        projection = projectToUnitSphere(center);
        totalError = totalError + norm(projection - center);
    end    
    avgCenterError = totalError / double(mesh.numF);
end

function avgEdgeLength = CalcAvgEdgeLength(mesh)
    [ii,jj,~] = find(mesh.Adjacency + mesh.Adjacency');
    edgeLengths = vecnorm(mesh.Vertices(ii,:) - mesh.Vertices(jj,:), 2, 2);
    avgEdgeLength = mean(edgeLengths);
    
end
