classdef MeshHW2 < MeshBasic
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        TriangleAreas
        VertexAreas
        InterpolantFtoV
        InterpolantVtoF
        
        % normals, for HW3
        FaceNormals
        VertexNormals
        EdgeNormalsMatrix % the E matrix
    end
    
    methods
        function obj = MeshHW2(filename)
           obj =  obj@MeshBasic(filename);
           obj = ComputeAreasAndInterpolant(obj);
           obj = ComputeNormals(obj);
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
        
        function obj = ComputeNormals(obj)
            % this is actually for HW3, but we put it here because it does
            % not rely on any differential operators, and it is fit for
            % everything
            
            % compute face normals
            i1 = obj.Faces(:,1);
            i2 = obj.Faces(:,2);
            i3 = obj.Faces(:,3);
            v1 = obj.Vertices(i1,:);
            v2 = obj.Vertices(i2,:);
            v3 = obj.Vertices(i3,:);
            obj.FaceNormals = cross(v1 - v2,v2 - v3,2);
            
            % compute vertex normals
            normals = obj.FaceNormals;
            normals = [interpolateFtoV(obj, normals(:,1)), interpolateFtoV(obj, normals(:,2)), interpolateFtoV(obj, normals(:,3)) ];
            normals = normals ./ vecnorm(normals, 2, 2);
            obj.VertexNormals = normals;
            
            % The E matrix
            
            % compute edge normals: vectors that are perpendicular to both
            % face normals and the edges
            n1 = cross(v2 - v3, obj.FaceNormals, 2);
            n2 = cross(v3 - v1, obj.FaceNormals, 2);
            n3 = cross(v1 - v2, obj.FaceNormals, 2);
            
            % make edge normals have same norm as their edges
            n1 = n1 ./ vecnorm(n1, 2, 2);
            n2 = n2 ./ vecnorm(n2, 2, 2);
            n3 = n3 ./ vecnorm(n3, 2, 2);
            
            n1 = n1 .* vecnorm(v2 - v3, 2, 2);
            n2 = n2 .* vecnorm(v3 - v1, 2, 2);
            n3 = n3 .* vecnorm(v1 - v2, 2, 2);

            % put them in the sparse matrix
            ii = repelem(obj.Faces, 1,3);
            jj = reshape(1:3*obj.numF, [obj.numF, 3]);
            jj = repmat(jj, [1,3]);
            vv = [n1 n2 n3];
            
            ii = reshape(ii, [3*3*obj.numF,1]);
            jj = reshape(jj, [3*3*obj.numF,1]);
            vv = reshape(vv, [3*3*obj.numF,1]);
            
            obj.EdgeNormalsMatrix = sparse(jj,ii,vv, 3*obj.numF,obj.numV);
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


