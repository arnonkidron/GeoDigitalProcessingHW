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
        
        % ComputeAreasAndInterpolant
        
        function obj = ComputeAreasAndInterpolant(obj)
            % compute triangle areas
            obj.TriangleAreas = zeros(obj.numF, 1);

            v1 = obj.Vertices(obj.Faces(:,1), :);
            v2 = obj.Vertices(obj.Faces(:,2), :);
            v3 = obj.Vertices(obj.Faces(:,3), :);
            obj.TriangleAreas = 0.5 * vecnorm(cross(v2-v1,v3-v1), 2, 2);
            
            % create VertexFaceAdjacency matrix, used to compute vertex areas
            ii = [1:obj.numF 1:obj.numF 1:obj.numF];
            jj = reshape(obj.Faces, [3 * obj.numF, 1]);
            VertexFaceAdjacency = sparse(ii,jj,1, obj.numF, obj.numV);
            
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
           obj.InterpolantVtoF = VertexFaceAdjacency / 3;
           
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
        
        %interpolantors
        function vertexFunction = interpolateFtoV(obj, faceFunction)
            product = obj.InterpolantFtoV * diag(faceFunction);
            vertexFunction = sum(product, 2);
        end
        
        function faceFunction = interpolateVtoF(obj, vertexFunction)
            product = obj.InterpolantVtoF * diag(vertexFunction);
            faceFunction = sum(product, 2);
        end
        
        
        % HW3: Differential Operators 
        
        function gradient = Gradient(obj, vertexFunc)
            gradient = zeros(obj.numF,3);
            for f = 1:obj.numF
                area = obj.TriangleAreas(f);
                i1 = obj.Faces(f,1);
                i2 = obj.Faces(f,2);
                i3 = obj.Faces(f,3);
                x1 = obj.Vertices(i1,:);
                x2 = obj.Vertices(i2,:);
                x3 = obj.Vertices(i3,:);
                p1 = Mesh.getProjection(x1, x2, x3);
                p2 = Mesh.getProjection(x2, x3, x1);
                p3 = Mesh.getProjection(x3, x1, x2);
                
         %       nablaB1 = (x1 - p1) / (2*area);
                nablaB2 = (x2 - p2) / (2*area);
                nablaB3 = (x3 - p3) / (2*area);
                
                % strange
                nablaB2 = nablaB2 / norm(nablaB2);
                nablaB2 = nablaB2 * norm(x1 - x3);
                
                nablaB3 = nablaB3 / norm(nablaB3);
                nablaB3 = nablaB3 * norm(x1 - x2);
                
                gradient(f,:) = ...
                    (vertexFunc(i2) - vertexFunc(i1)) * nablaB2 + ...
                    (vertexFunc(i3) - vertexFunc(i1)) * nablaB3;
            end
        end
        
        
    end
    methods(Static)
        
        function p = getProjection(x1, x2, x3)
            % input: 3 vertices of a triangle, counterclockwise
            % output: the projection of x1 unto x2x3
            u = x1 - x2;
            v = x3 - x2;
            p = x2 + dot(u, v) / dot(v, v) * v;
        end   
    end
end


