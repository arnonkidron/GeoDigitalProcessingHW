classdef MeshHW2 < MeshBasic
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        TriangleAreas
        VertexAreas
        InterpolantFtoV
        InterpolantVtoF
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
            % compute ...
            TriangleAreasRepeatedInverse = power(obj.TriangleAreas,-1);
            TriangleAreasRepeatedInverse = repmat(TriangleAreasRepeatedInverse, 3, 1);
            TriangleAreasRepeatedInverse = diag(TriangleAreasRepeatedInverse);
            
            % compute the E matrix, EdgeNormals
            ii = [1:obj.numF;(obj.numF+1):2*obj.numF;(2*obj.numF+1):3*obj.numF];
            ii = repelem(ii, 1,3);
            ii = reshape(ii, [3*3*obj.numF, 1]);
            
            jj = repelem(obj.Faces, 1,3);
%            jj = [obj.Faces(:,1) obj.Faces(:,1) obj.Faces(:,1), obj.Faces(:,2) obj.Faces(:,2) obj.Faces(:,2),obj.Faces(:,3) obj.Faces(:,3) obj.Faces(:,3)];
            jj = reshape(jj', [3*3*obj.numF,1]);
            
            vv = zeros(3*3*obj.numF,1);
            
%             EdgeNormals = sparse(3*obj.numF, obj.numV);
            % TODO: vectorize
            index = 0;
            for f = 1:obj.numF
%                 i1 = obj.Faces(f,1);
%                 i2 = obj.Faces(f,2);
%                 i3 = obj.Faces(f,3);
                 
                x1 = obj.Vertices(obj.Faces(f,1),:);
                x2 = obj.Vertices(obj.Faces(f,2),:);
                x3 = obj.Vertices(obj.Faces(f,3),:);
                p1 = MeshHW2.getProjection(x1, x2, x3);
                p2 = MeshHW2.getProjection(x2, x3, x1);
                p3 = MeshHW2.getProjection(x3, x1, x2);
                
                Je1 = x1 - p1;
                Je2 = x2 - p2;
                Je3 = x3 - p3;
                
                Je1 = Je1 / norm(Je1);
                Je2 = Je2 / norm(Je2);
                Je3 = Je3 / norm(Je3);
                
                Je1 = Je1 * norm(x3-x2);
                Je2 = Je2 * norm(x3-x1);
                Je3 = Je3 * norm(x1-x2);
                
                vv(index + 1:index + 9) = [Je1 Je2 Je3];
                index = index + 9;
                
%                 EdgeNormals(f,              i1) = Je1(1);
%                 EdgeNormals(f + obj.numF,   i1) = Je1(2);
%                 EdgeNormals(f + 2*obj.numF, i1) = Je1(3);
%                 
%                 EdgeNormals(f,              i2) = Je2(1);
%                 EdgeNormals(f + obj.numF,   i2) = Je2(2);
%                 EdgeNormals(f + 2*obj.numF, i2) = Je2(3);
%                 
%                 EdgeNormals(f,              i3) = Je3(1);
%                 EdgeNormals(f + obj.numF,   i3) = Je3(2);
%                 EdgeNormals(f + 2*obj.numF, i3) = Je3(3);
            end
            
%             fast = sparse(ii,jj,vv, 3*obj.numF, obj.numV);
            
            EdgeNormals = sparse(ii,jj,vv, 3*obj.numF, obj.numV);
            
            grad = sparse(TriangleAreasRepeatedInverse * EdgeNormals) / 2;
            
            % transpose if necessary
            if(size(vertexFunc,1) == 1)
                vertexFunc = vertexFunc';
            end
            
            gradient = grad * vertexFunc;
        end
        
        function gradient = OldGradient(obj, vertexFunc)
            gradient = zeros(obj.numF,3);
            for f = 1:obj.numF
                area = obj.TriangleAreas(f);
                i1 = obj.Faces(f,1);
                i2 = obj.Faces(f,2);
                i3 = obj.Faces(f,3);
                x1 = obj.Vertices(i1,:);
                x2 = obj.Vertices(i2,:);
                x3 = obj.Vertices(i3,:);
                p1 = MeshHW2.getProjection(x1, x2, x3);
                p2 = MeshHW2.getProjection(x2, x3, x1);
                p3 = MeshHW2.getProjection(x3, x1, x2);
                
                nablaB1 = (x1 - p1) / (2*area);
                nablaB2 = (x2 - p2) / (2*area);
                nablaB3 = (x3 - p3) / (2*area);
                
                % same for nablaB1
                nablaB1 = nablaB1 / norm(nablaB1);
                nablaB1 = nablaB1 * norm(x3 - x2);
                
                % make nablaB2 the same size as the opposite edge
                nablaB2 = nablaB2 / norm(nablaB2);
                nablaB2 = nablaB2 * norm(x1 - x3);
                
                % same for nablaB3
                nablaB3 = nablaB3 / norm(nablaB3);
                nablaB3 = nablaB3 * norm(x1 - x2);

                gradient(f,:) = ...
                    vertexFunc(i1) * nablaB1 + ...
                    vertexFunc(i2) * nablaB2 + ...
                    vertexFunc(i3) * nablaB3;
                
            end
        end
        
        function div = Divergence(obj, vectorField)
            GF = diag(obj.TriangleAreas); % replicate
            GV = diag(obj.VertexAreas);
            grad = Gradient(obj, vectorField);
            
            
            div = GF * grad * GV;
            
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


