classdef Mesh < MeshHW2
    % Mesh for HW3, with differential geometry stuff
    
    properties
        Grad % 3|F| x |V|
        Div  % |V| x 3|F|
        Laplacian  % |V| x |V|, computed as minus Divergence of the Gradient
        CotLaplacian  % computed by the cotangent formula
        
        MeanCurvature
        GaussianCurvature
        
        CotangentWeights % for analysis1
    end
    
    methods
        function obj = Mesh(filename)
            obj =  obj@MeshHW2(filename);
            obj = ComputeDifferentialOperators(obj);
            obj = ComputeCurvatures(obj);
        end
        
        function obj = ComputeDifferentialOperators(obj)
            % compute GF & GV, the area matrices
            % Note that GF is repeated 3 times
            TriangleAreasRepeatedInverse = power(obj.TriangleAreas,-1);
            TriangleAreasRepeatedInverse = repmat(TriangleAreasRepeatedInverse, 3, 1);
            TriangleAreasRepeatedInverse = sparse(diag(TriangleAreasRepeatedInverse));
            
            TriangleAreasRepeated = repmat(obj.TriangleAreas, 3, 1);
            TriangleAreasRepeated = sparse(diag(TriangleAreasRepeated));
            
            VertexAreasInverse = power(obj.VertexAreas, -1);
            VertexAreasInverse = sparse(diag(VertexAreasInverse));
            
            % Gradient
            EdgeNormals = ComputeEdgeNormals(obj);
%             EdgeNormals = obj.EdgeNormalsMatrix;
            obj.Grad = sparse(TriangleAreasRepeatedInverse * EdgeNormals) / 2;
            
            % Divergence
            obj.Div = -VertexAreasInverse * obj.Grad' * TriangleAreasRepeated;
           
            % Laplacian
            obj.Laplacian = -obj.Div * obj.Grad;
            obj.CotLaplacian = LaplacianByCotFormula(obj);
            
            % cotangent weights
            obj.CotangentWeights = ...
                EdgeNormals' * ...
                TriangleAreasRepeatedInverse * ...
                EdgeNormals/ 4;
        end
        
        function EdgeNormals = ComputeEdgeNormals(obj)
            ii = [1:obj.numF;(obj.numF+1):2*obj.numF;(2*obj.numF+1):3*obj.numF];
            ii = repelem(ii, 1,3);
            ii = reshape(ii, [3*3*obj.numF, 1]);
            
            jj = repelem(obj.Faces, 1,3);
            jj = reshape(jj', [3*3*obj.numF,1]);
            
            vv = zeros(3*3*obj.numF,1);
            
            % TODO: vectorize
            index = 0;
            for f = 1:obj.numF
                x1 = obj.Vertices(obj.Faces(f,1),:);
                x2 = obj.Vertices(obj.Faces(f,2),:);
                x3 = obj.Vertices(obj.Faces(f,3),:);
                p1 = Mesh.getProjection(x1, x2, x3);
                p2 = Mesh.getProjection(x2, x3, x1);
                p3 = Mesh.getProjection(x3, x1, x2);
                
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
               
                % % what we really meant to do
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
            
            EdgeNormals = sparse(ii,jj,vv, 3*obj.numF, obj.numV);
        end
        
        function obj = ComputeCurvatures(obj)
            % Mean Curvature - very simple
            obj.MeanCurvature = vecnorm(obj.Laplacian, 2,2);
            
            % Gaussian Curvature
            
            % put the angles in a |F| x |V| matrix
            A = ComputeAngles(obj);
            ii = repmat([1:obj.numF], [1,3]);
            jj = obj.Faces;
            B = sparse(ii, jj, A, obj.numF, obj.numV);
            
            % compute 
            sumAngles = sum(B, 1)';
            curv = 2*pi - sumAngles;
            obj.GaussianCurvature = curv ./ obj.VertexAreas;
        end
        
        function Gradient = CalcGradient(obj, vertexFunc)
            % transpose if necessary
            if(size(vertexFunc,1) == 1)
                vertexFunc = vertexFunc';
            end
            
            Gradient = obj.Grad * vertexFunc;
        end
        
        function Gradient = OldGradient(obj, vertexFunc)
            Gradient = zeros(obj.numF,3);
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

                Gradient(f,:) = ...
                    vertexFunc(i1) * nablaB1 + ...
                    vertexFunc(i2) * nablaB2 + ...
                    vertexFunc(i3) * nablaB3;
                
            end
        end
        
        function Divergence = CalcDivergence(obj, faceVectorField)
            if(size(faceVectorField) == [obj.numF, 3])
                faceVectorField = reshape(faceVectorField, [3*obj.numF, 1]);
            end
            Divergence = obj.Div * faceVectorField;
        end
        
        function lap = CalcLaplacian(obj, vertexFunc)
            lap = obj.Laplacian * vertexFunc;            
        end
        
        function A = ComputeAngles(obj)
            % @returns the matrix A contains all angles
            % the 1st column contains the angle of the 1st vertex in every face, etc.
            
            % compute all edge lengths
            v1 = obj.Vertices(obj.Faces(:,1),:);
            v2 = obj.Vertices(obj.Faces(:,2),:);
            v3 = obj.Vertices(obj.Faces(:,3),:);
            L1 = vecnorm(v2 - v3, 2, 2);
            L2 = vecnorm(v1 - v3, 2, 2);
            L3 = vecnorm(v2 - v1, 2, 2);
            
            % compute cosines of all angles, using cosine theorem
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
            
            A = [A1, A2, A3];
            
            % finally, apply arccos
            A = acos(A);
        end
        
        function laplaceMatrix = LaplacianByCotFormula(obj)
            A = ComputeAngles(obj);
            
            % cotangent formula
            I = [obj.Faces(:,1); obj.Faces(:,2); obj.Faces(:,3)];
            J = [obj.Faces(:,2); obj.Faces(:,3); obj.Faces(:,1)];
            S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);
            In = [I;J;I;J];
            Jn = [J;I;I;J];
            Sn = [-S;-S;S;S];
            
            laplaceMatrix = sparse(In,Jn,Sn, obj.numV,obj.numV); 
        end
        
        function [fig, p] = RenderGradient(obj, vertexFunc)
            g = CalcGradient(obj, vertexFunc);
            [fig, p] = RenderVectorField(obj, vertexFunc, g);
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