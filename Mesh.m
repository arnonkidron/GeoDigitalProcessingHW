classdef Mesh
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        Vertices
        Faces
        Adjacency % as directed edges
        numV
        numE
        numF
        TriangleAreas
        VertexAreas
        InterpolantFtoV
        InterpolantVtoF
    end
    
    methods
        % TODO: Mesh([]) opens a dialog to choose file
        function obj = Mesh(filename)
            % input: path for an OFF file to read from
            % output: a Mesh object 
            
            %open file
            file = fopen(filename);
            
            % skip the first line that says "OFF"
            fgetl(file);
            
            % read the number of vertices, faces into numV, numF
            tmpCell = deal(textscan(file, "%d %d %d", 1));
            [obj.numV, obj.numF, ~] = deal(tmpCell{:});
            
            % read the vertices into the matrix 
            tmpCell = textscan(file, "%f %f %f", obj.numV);
            [a,b,c] = deal(tmpCell{:});
            obj.Vertices = [a b c];
            
            % read the faces
            tmpCell = textscan(file, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d", obj.numF, ...
            'CollectOutput',   true, 'EmptyValue', -1 );
            tmpMat = tmpCell{:};
            sizesCol = tmpMat(:,1);
            maxFaceSize = max(sizesCol);
            minFaceSize = min(sizesCol);
            tmpMap = tmpMat(:,2:1+maxFaceSize);
            
            % shift indices by 1 to fit matlab indexing
            tmpMap = tmpMap + 1;
            obj.Faces = tmpMap;            

            % close file
            fclose(file);
            
            % create adjacency matrix
            if(minFaceSize < maxFaceSize)
                disp("Too hard")
                return
            end
            v1 = reshape(obj.Faces, [obj.numF * maxFaceSize, 1]);
            v2 = circshift(v1, -obj.numF);
            
            obj.Adjacency = sparse(v1, v2, 1, obj.numV,obj.numV);
            obj.numE = nnz(obj.Adjacency + obj.Adjacency') / 2;
            
            % more properties
            obj = ComputeAreasAndInterpolant(obj);
            
        end
        
        function Write(obj, outputFilename)
            % Exports the Mesh to an OFF file 
            fid = fopen(outputFilename, 'wt');
            fprintf(fid, "OFF\n");
            fprintf(fid, "%d %d %d\n", obj.numV, obj.numF, 0);
            fprintf(fid, "%g %g %g\n", obj.Vertices');
            fprintf(fid, "%d %d %d %d\n", [(3 * ones(obj.numF, 1))  (obj.Faces - 1)]');
            fclose(fid);
            
            
        end
        
        function Render(obj, colors)
            % opens a figure of the mesh
            
            if(isempty(colors))
                colors = "Green";
            end
            
            figure
            colorbar
            p = patch('Faces', obj.Faces, 'Vertices', obj.Vertices);
            
            numC = size(colors,1);
            if(numC == obj.numV) % colors for each vertex
                p.FaceVertexCData = colors;
                p.FaceColor = 'interp';
            elseif(numC == obj.numF) % colors for each face
                p.FaceVertexCData = colors;
                p.FaceColor = 'flat';
            else % single color for whole mesh
                p.FaceColor = colors;
            end
            view(3);
            
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
        
        % topology measures
        
        function genus = CalcGenus(obj)
            b = CalcBoundaryComponents(obj);
            chi = obj.numV - obj.numE + obj.numF;
            genus = (2 - chi - b)/2;
        end
        
        function boundaryEdges = CalcBoundaryEdges(obj)
            diff_matrix = obj.Adjacency - transpose(obj.Adjacency);
            boundaryEdges = nnz(diff_matrix)/2;
        end
        
        function boundaryCC = CalcBoundaryComponents(obj)
            diff_matrix = obj.Adjacency - transpose(obj.Adjacency);
            diff_matrix = abs(diff_matrix);
            G = graph(diff_matrix);
            [~, binsizes] = conncomp(G);
            boundaryCC = nnz(binsizes > 1);
        end
        
    end
end


