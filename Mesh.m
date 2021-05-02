classdef Mesh < MeshWithoutAreaStuff
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        TriangleAreas
        VertexAreas
        InterpolantFtoV
        InterpolantVtoF
    end
    
    methods
        function obj = Mesh(filename)
           obj =  obj@MeshWithoutAreaStuff(filename);
           obj = ComputeAreasAndInterpolant(obj);
        end
        
        function obj = ComputeAreasAndInterpolant(obj)
            % compute triangle areas
            obj.TriangleAreas = zeros(obj.numF, 1);
            
            for f = 1:obj.numF
                v1 = obj.Vertices(obj.Faces(f,1), :);
                v2 = obj.Vertices(obj.Faces(f,2), :);
                v3 = obj.Vertices(obj.Faces(f,3), :);
                
                obj.TriangleAreas(f) = 0.5 * norm(cross(v2-v1,v3-v1));                
            end
            
            % create VertexFaceAdjacency matrix, used to compute vertex areas
            ii = [1:obj.numF 1:obj.numF 1:obj.numF];
            jj = reshape(obj.Faces, [3 * obj.numF, 1]);
            VertexFaceAdjacency = sparse(ii,jj,1);
            
            % compute vertex areas
            
            product = diag(obj.TriangleAreas) * VertexFaceAdjacency;
            obj.VertexAreas = sum(product, 1)' / 3;  
           
           % compute interplant           
           obj.InterpolantFtoV = ...
               sparse(diag(obj.VertexAreas .^ -1)) * ...
               VertexFaceAdjacency' * ...
               sparse(diag(obj.TriangleAreas));
           obj.InterpolantFtoV = obj.InterpolantFtoV / 3;
           
           % compute other interpolant
           obj.InterpolantVtoF = ...
               sparse(diag(obj.TriangleAreas .^ -1)) * ...
               obj.InterpolantFtoV' * ...
               sparse(diag(obj.VertexAreas));
            
        end
        
        function triangleAreas = GetTriangleAreas(obj)
            triangleAreas = obj.TriangleAreas;
        end
        
        function vertexAreas = GetVertexAreas(obj)
            vertexAreas = obj.VertexAreas;
        end
        
        %interpolantors
        function vertexFunction = interpolateFtoV(obj, faceFunction)
            product = obj.InterpolantFtoV * diag(faceFunction);
            vertexFunction = sum(product, 2);
        end
        
        function faceFunction = interpolateVtoF(obj, vertexFunction)
            product = obj.InterpolantVtoF * diag(vertexFunction);
            faceFunction = sum(product, 2);
        end
        
    end
end


