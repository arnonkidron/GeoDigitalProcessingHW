classdef MeshHW2 < MeshBasic
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        TriangleAreas
        VertexAreas
        InterpolantFtoV
        InterpolantVtoF
        VertexFaceAdjacency
    end
    
    methods
        function obj = MeshHW2(filename)
           obj =  obj@MeshBasic(filename);
           obj = ComputeAreasAndInterpolant(obj);
        end
        
        % ComputeAreasAndInterpolant
        
        function obj = ComputeAreasAndInterpolant(obj)
            % compute triangle areas
            obj.TriangleAreas = zeros(obj.numF, 1);

            v1 = obj.Vertices(obj.Faces(:,1), :);
            v2 = obj.Vertices(obj.Faces(:,2), :);
            v3 = obj.Vertices(obj.Faces(:,3), :);
            obj.TriangleAreas = 0.5 * vecnorm(cross(v2-v1,v3-v1), 2, 2);
            
            % create obj.VertexFaceAdjacency matrix, used to compute vertex areas
            ii = [1:obj.numF 1:obj.numF 1:obj.numF];
            jj = reshape(obj.Faces, [3 * obj.numF, 1]);
            obj.VertexFaceAdjacency = sparse(ii,jj,1, obj.numF, obj.numV);
            
            % compute vertex areas
            
            product = diag(obj.TriangleAreas) * obj.VertexFaceAdjacency;
            obj.VertexAreas = sum(product, 1)' / 3;  
           
           % compute interplant           
           obj.InterpolantFtoV = ...
               sparse(diag(obj.VertexAreas .^ -1)) * ...
               obj.VertexFaceAdjacency' * ...
               sparse(diag(obj.TriangleAreas));
           obj.InterpolantFtoV = obj.InterpolantFtoV / 3;
           
           % compute other interpolant
           obj.InterpolantVtoF = obj.VertexFaceAdjacency / 3;
           
           % naive calculation by the HW2 instructions
%            obj.InterpolantVtoF = ...
%                sparse(diag(obj.TriangleAreas .^ -1)) * ...
%                obj.InterpolantFtoV' * ...
%                sparse(diag(obj.VertexAreas));
           
           
            
        end
        
        function triangleAreas = GetTriangleAreas(obj)
            triangleAreas = obj.TriangleAreas;
        end
        
        function vertexAreas = GetVertexAreas(obj)
            vertexAreas = obj.VertexAreas;
        end
        
        %interpolators
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


