classdef MeshPCA
    properties
        Mesh
        Matrix
        Eigenvalues
        Eigenvectors
    end
    
    methods
        function obj = MeshPCA(url, max_k)
            obj.Mesh = Mesh(url);
            obj.Matrix = obj.Mesh.CotangentWeights;
            
            [obj.Eigenvectors, obj.Eigenvalues, failure] = ...
                eigs(obj.Matrix, max_k, 'lm');
            if(failure)
                disp("Warning: not all eigenvalues have converged");
            end
        end
        
        function fig = Render9Eigenfunctions(obj)
            fig = RenderSeveralFunctions(obj.Mesh, obj.Eigenvectors(:,1:9), []);
        end
        
        function f = GetHatFunction(obj, i)
            f = zeros(obj.Mesh.numV, 1);
            f(i) = 1;
        end
        
        function f = GetMaximalEigenfunction(obj)
            [f, ~] = eigs(obj.Matrix, 1, 'sm');
        end
        
        function [g, error] = ProjectAndBack(obj, f, k)
            B = obj.Eigenvectors(:,1:k);
            g = B * (B' * f);
            error = norm(g - f);
        end
    end
end