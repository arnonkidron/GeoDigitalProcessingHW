classdef Mesh
    % a 3D Mesh
    % implemented with a Shared Vertex Table & Vertex Adjacency Matrix
    
    properties
        Vertices
        Faces
        Adjacency
        TriangleAreasResult
    end
    
    methods
        function obj = Mesh(filename)
            % input: path for an OFF file to read from
            % output: a Mesh object 
            
            obj.TriangleAreasResult = [];
            
            %open file
            file = fopen(filename);
            
            % skip the first line that says "OFF"
            fgetl(file);
            
            
            
            % read the number of vertices, faces & edges
            tmpCell = deal(textscan(file, "%d %d %d", 1));
            [numV, numF, ~] = deal(tmpCell{:});
            
            % read the vertices into the matrix 
            tmpCell = textscan(file, "%f %f %f", numV);
            [a,b,c] = deal(tmpCell{:});
            obj.Vertices = [a b c];
            
            % read the faces
            tmpCell = textscan(file, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d", numF, ...
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
            v1 = reshape(obj.Faces, [numF * maxFaceSize, 1]);
            v2 = circshift(v1, -numF);
            
            obj.Adjacency = sparse(v1, v2, 1, numV,numV);
            obj.Adjacency = obj.Adjacency + transpose(obj.Adjacency);
        end
        
        function outputArg = Write(obj, outputFilename)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
        end
        
        function Render(obj, faceColor)
            % opens a figure of the mesh
            
            if(isempty(faceColor))
                faceColor = "Green";
            end
            
%             c = obj.Vertices(:,1);
%             c = c - min(c);
%             c = c / max(c);
%             c = [c c c];
            
            patch('Faces', obj.Faces, 'Vertices', obj.Vertices, ...
                'FaceColor', faceColor);
            view(3);
            
        end
        
        function triangleAreas = TriangleAreas(obj)
            if(~isempty(obj.TriangleAreasResult))
                triangleAreas = obj.TriangleAreasResult;
                return
            end
            
            numF = size(obj.Faces, 1);
            obj.TriangleAreasResult = zeros(numF, 1);
            
            for f = 1:numF
                v1 = obj.Vertices(obj.Faces(f,1), :);
                v2 = obj.Vertices(obj.Faces(f,2), :);
                v3 = obj.Vertices(obj.Faces(f,3), :);
                
                obj.TriangleAreasResult(f) = 0.5 * norm(cross(v2-v1,v3-v1));                
            end            
            
            triangleAreas = obj.TriangleAreasResult;
        end
        
        function vertexAreas = VertexAreas(obj)
            triangleAreas = TriangleAreas(obj);
            numF = size(obj.Faces, 1);
            numV = size(obj.Vertices, 1);
            
            VertexFaceAdjacency = zeros(numF, numV);
            for v = 1:numV
                VertexFaceAdjacency(:,v) = any(obj.Faces == v, 2);
            end
            VertexFaceAdjacency = diag(triangleAreas) * ...
                VertexFaceAdjacency;
            
            vertexAreas = sum(VertexFaceAdjacency, 2);
        end
        
        % private functions
%         function computeEdgeLengths(obj)
%             for face = 1:size(obj.Faces,1)
%                 for v = 1:3
%                     u = v + 1;
%                     if u == 4
%                         u = 1;
%                     end
%                     vIndex = obj.Faces(f,v);
%                     uIndex = obj.Faces(f,u);
%                     e = obj.Vertices(vIndex, :) - obj.Vertices(uIndex, :);
%                     obj.FacesEdgeLengths(v, face) = norm(e);
%                 end
%             end
%         end
    end
end


