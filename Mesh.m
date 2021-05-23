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
            obj.Grad = sparse(TriangleAreasRepeatedInverse * obj.EdgeNormalsMatrix) ...
                / 2;
            
            % Divergence
            obj.Div = -VertexAreasInverse * obj.Grad' * TriangleAreasRepeated;
           
            % Laplacian
            obj.Laplacian = -obj.Div * obj.Grad;
            obj.CotLaplacian = LaplacianByCotFormula(obj);
            
            % cotangent weights
            obj.CotangentWeights = ...
                obj.EdgeNormalsMatrix' * ...
                TriangleAreasRepeatedInverse * ...
                obj.EdgeNormalsMatrix ...
                / 4;
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
        
        function [fig, p, arrows] = RenderGradient(obj, vertexFunc)
            g = CalcGradient(obj, vertexFunc);
            [fig, p, arrows] = RenderVectorField(obj, vertexFunc, g);
        end
        
        
    end
end