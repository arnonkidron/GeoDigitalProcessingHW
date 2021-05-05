classdef MeshWithoutAreaStuff
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        Name
        Vertices
        Faces
        Adjacency % as directed edges
        numV
        numE
        numF
    end
    
    methods
        function obj = MeshWithoutAreaStuff(filename)
            % input: path for an OFF file to read from
            % output: a Mesh object 
            
            %open file
            if(isempty(filename))
                filename = askUserForSingleOFFfile();
            end
            file = fopen(filename);
            
            [~, obj.Name, ~] = fileparts(filename);
            
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
        
        function [fig, p] = Render(obj, colors)    
            % @param colors: can be either one of the following
            % * a single color or [] - to paint all faces with that color
            % * a numeric vector of size obj.numF - to paint each face by its value
            % * a numeric vector of size obj.numV - to paint each vertex by its
            % value, and interpolate the fragment colors from the vertices
        
            % opens a figure of the mesh
            
            if(isempty(colors))
                colors = "Green";
            end
            
            fig = figure;
            set(fig, 'Name', obj.Name);
            colorbar;
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
        
        function [fig, p] = RenderWireframe(obj, vertexColors)
            % opens a figure of the mesh
            
            if(isempty(vertexColors))
                vertexColors = "Green";
            end
            
            fig = figure;
            set(fig, 'Name', obj.Name);
            colorbar;
            p = patch('Faces', obj.Faces, 'Vertices', obj.Vertices);
            colormap parula;
            
            p.FaceColor = 'c';
            p.FaceAlpha = 0.1;
            
            
            numC = size(vertexColors,1);
            if(numC == obj.numV) % colors for each vertex
                p.FaceVertexCData = vertexColors;
                p.EdgeColor = 'interp';
            else % just one color
                p.EdgeColor = vertexColors;
            end
            
            view(3);
            
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
        
        
        function G = GetWeightedGraph(obj)
            % returns the graph of vertices, where two vertices are
            % neighbours iff they have a common edge
            % the weight are the vertex distances
            [ii,jj,~] = find(obj.Adjacency + obj.Adjacency');
            ww = vecnorm(obj.Vertices(ii,:) - obj.Vertices(jj,:), 2, 2);
            G = graph(sparse(ii,jj,ww, obj.numV,obj.numV), 'lower');
            
        end
    end
end


