mesh = Mesh("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\disk.off");
vertexFunc = GetValences(mesh);
g = Gradient(mesh, vertexFunc);
RenderVectorField(mesh, vertexFunc, g);

view(2)
colormap flag

