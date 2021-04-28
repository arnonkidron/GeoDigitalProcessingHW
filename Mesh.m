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
        TriangleAreasResult
    end
    
    methods
        % TODO: Mesh([]) opens a dialog to choose file
        function obj = Mesh(filename)
            % input: path for an OFF file to read from
            % output: a Mesh object 
            
            obj.TriangleAreasResult = [];
            
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
        
        function triangleAreas = TriangleAreas(obj)
            % returns the vector of triangle areas
            if(~isempty(obj.TriangleAreasResult))
                triangleAreas = obj.TriangleAreasResult;
                return
            end
            
            obj.TriangleAreasResult = zeros(obj.numF, 1);
            
            for f = 1:obj.numF
                v1 = obj.Vertices(obj.Faces(f,1), :);
                v2 = obj.Vertices(obj.Faces(f,2), :);
                v3 = obj.Vertices(obj.Faces(f,3), :);
                
                obj.TriangleAreasResult(f) = 0.5 * norm(cross(v2-v1,v3-v1));                
            end            
            
            triangleAreas = obj.TriangleAreasResult;
        end
        
        function vertexAreas = VertexAreas(obj)
            % return the vector of vertex areas, defined as 1/3 of their
            % adjacent triangles' areas
            triangleAreas = TriangleAreas(obj);
            
            VertexFaceAdjacency = zeros(obj.numF, obj.numV);
            for v = 1:obj.numV
                VertexFaceAdjacency(:,v) = any(obj.Faces == v, 2);
            end
            VertexFaceAdjacency = diag(triangleAreas) * ...
                VertexFaceAdjacency;
            
            vertexAreas = sum(VertexFaceAdjacency, 1)' / 3;
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
        
        %interpolants
        function fFunc = interpolateVtoF(obj, vFunc)
            
        end
        
        function vFunc = interpolateFtoV(obj, fFunc)
            
        end
        
        
        
    end
end


