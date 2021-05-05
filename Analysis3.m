urls = askUserForMultipleOFFfiles();

for url_1 = urls
    url = url_1{1};
    mesh = Mesh(url);
    nullSpace = null(full(mesh.InterpolantFtoV));
    if(isempty(nullSpace))
        [~,name,~] = fileparts(url);
        disp("Zero kernel for " + name);
        continue;
    end
    if(mesh.numV >= 1000)
        disp("Skipping " + name + " because it has too many vertices");
    end
    
    limit = 10;
    i = 0;
    for vertexFunction = nullSpace
        i = i + 1;
        if i >= limit
            break
        end
        faceFunction = interpolateFtoV(mesh, vertexFunction);
        if(min(faceFunction) == 0 && max(faceFunction) == 0)
            disp("Failure")
        end
        
        fig = Render(mesh, faceFunction);
        
        
        
        % adjust colorbar to have 0 in the middle
        cb = colorbar;
        oldLimits = get(cb, 'Limits');
        l = max(abs(oldLimits));
        set(cb, 'Limits', [-l l]);
        colormap hot
    end
    if(~isempty(nullSpace))
        break
    end
end



