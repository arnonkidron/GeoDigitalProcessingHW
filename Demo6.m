
dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
filename = "disk.off";
mesh = Mesh(dir + filename);

hat = zeros(mesh.numV,1);
hat(226) = 1;
myRender(mesh, hat, "a hat function on vertices");

fHat = interpolateVtoF(mesh, hat);
myRender(mesh, fHat, "an interpolated function for faces");

reconstructedHat = interpolateFtoV(mesh, fHat);
myRender(mesh, reconstructedHat, "interpolated back");

function myRender(mesh, f, titleText)
    fig = Render(mesh, f);
    view(2);
    set(get(fig, 'CurrentAxes'), 'Clim', [0 1]);
    cb = colorbar;
    set(cb, 'ticks', unique(f));
    set(cb, 'ylim', [min(f) max(f)]);
    colormap jet
    
    title(titleText);
end