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
    for faceFunction = nullSpace
        i = i + 1;
        if i >= limit
            break
        end
        
        % normalize
        mag = max(abs(faceFunction));
        faceFunctionNormalized = faceFunction / mag;
        
        % render
        fig = Render(mesh, faceFunctionNormalized);
        view(2);
        
        
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



