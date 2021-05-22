classdef MeshBasic
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
        function obj = MeshBasic(filename)
            % input: path for an OFF file to read from
            % output: a MeshHW2 object
            
            %open file
            if(isempty(filename))
                filename = askUserForSingleOFFfile();
            end
            file = fopen(filename);
            if(file == -1)
                error("Fail to open url");
            end
            
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
            % Exports the MeshHW2 to an OFF file
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
            % @returns: the figre, and the patch
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
                set(colorbar,'visible','off')
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
        
        function [fig, p, arrows] = RenderVectorField(obj, func, vectorField)
            % @returns: the figure, the patch
            % and the quiver object so that the user may change the arrows color as he wish
            [fig,p] = Render(obj, func);
          %  alpha clear
            hold on
            sizeVectorField = size(vectorField,1);
            
            if(sizeVectorField == 3*obj.numF || sizeVectorField == 3*obj.numV)
                vectorField = reshape(vectorField, [sizeVectorField / 3, 3]);
                sizeVectorField = size(vectorField,1);
            end
            
            if any(abs(vecnorm(vectorField ,2,2) - 1) > 0.000001)
%                disp("wrong");
            end
            
            
            if(sizeVectorField == obj.numV)
                % function on vertices
                arrows = quiver3(...
                    obj.Vertices(:,1), ...
                    obj.Vertices(:,2), ...
                    obj.Vertices(:,3), ...
                    vectorField(:,1), ... 
                    vectorField(:,2), ... 
                    vectorField(:,3), ... 
                    'color', "Red" ...
                    );
            elseif(sizeVectorField == obj.numF)
                % function on faces
                centers = getTriangleCenters(obj);
                arrows = quiver3(...
                    centers(:,1), ...
                    centers(:,2), ...
                    centers(:,3), ...
                    vectorField(:,1), ... 
                    vectorField(:,2), ... 
                    vectorField(:,3), ... 
                    'color', "Red", ...
                    'linewidth', 1 ...
                    );
            else
                disp("Incsompatible vector field!")
            end    
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
        
        function valences = GetValences(obj)
            A = obj.Adjacency + obj.Adjacency';
            valences = full(sum(A ~= 0, 2));
        end
        
        function center = getTriangleCenter(obj, f)
            % input: a mesh
            %        a face index
            % output: its center point
            v1 = obj.Vertices(obj.Faces(f, 1), :);
            v2 = obj.Vertices(obj.Faces(f, 2), :);
            v3 = obj.Vertices(obj.Faces(f, 3), :);

            center = mean([v1; v2; v3], 1);
        end
        
        function centers = getTriangleCenters(obj)
            % input: a mesh
            % output: the center points of all faces
            centers = zeros(obj.numF, 3);
            for f = 1:obj.numF
                centers(f,:) = getTriangleCenter(obj, f);
            end
        end
    end
end


