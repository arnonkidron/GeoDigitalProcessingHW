MeshHW2 = MeshHW2("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\disk.off");
 
vertexFunc = GetValences(MeshHW2);
% RenderGradient(MeshHW2, vertexFunc);
% 
% view(2)
% colormap flag
% 



%MeshHW2 = MeshHW2("C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\triangle.off");
%vertexFunc = [0 0 1]';
simple = OldGradient(MeshHW2, vertexFunc);
complex = reshape(Gradient(MeshHW2, vertexFunc), [MeshHW2.numF, 3]);
RenderVectorField(MeshHW2, vertexFunc, simple);
view(2);
colormap flag;
RenderVectorField(MeshHW2, vertexFunc, complex);
view(2);
colormap flag;