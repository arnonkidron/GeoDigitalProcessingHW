dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
meshName = "sphere_s2.off";
global EDGE_ALPHA;
EDGE_ALPHA = 0.05;

mesh = MeshBasic(dir + meshName);
[~, patches] = RenderSphericalHarmonics(mesh);

% meshSmoother = MeshSmoother(dir + meshName, max(k));
% [~, patches] = RenderSomeEigenfunctions(meshSmoother, 25, EDGE_ALPHA);

function [fig, patches] = RenderSphericalHarmonics(mesh)
    global EDGE_ALPHA;
    sphereCoords = getSphericalCoordinates(mesh.Vertices);
    max_degree = 4;
    A = ones(mesh.numV, (max_degree+1) ^2);
    titles = cell((max_degree+1) ^2, 1);
    for l=0:max_degree
        for m=0:l
            theta = sphereCoords(:,1);
            phi = sphereCoords(:,2);
            harmonic = Ylm(l,m,theta,phi, false);

            index = m + l*(max_degree+1) + 1;
            A(:, index) = harmonic;
            
            titles{index} = sprintf("l=%d, m=%d", l, m);
        end
        for m=l+1:max_degree
            index = m + l*(max_degree+1) + 1;
            A(:,index) = NaN;
        end
    end
    
    [fig,patches] = RenderSeveralFunctions(mesh, A, titles, EDGE_ALPHA);
end

% function adjustEdgeAlpha(patches, val)
%     for i = 1:numel(patches)
%         p = patches{i};
%         if isempty(p)
%             continue;
%         end
%         
%         set(p, 'EdgeAlpha', val);
%     end
% end

function sphereCoords = getSphericalCoordinates(cartesianCoords)
    x = cartesianCoords(:,1);
    y = cartesianCoords(:,2);
    z = cartesianCoords(:,3);
    theta = atan(sqrt(x.^2 + y.^2) ./ z);
    phi = atan(y ./ x);
    phi(isnan(phi))=0;
    sphereCoords = [theta, phi];
end

function Ylm = Ylm(l,m,theta,phi, shouldNormalize)
    theta = theta(:);
    phi = phi(:)';
    
    % compute the legendre polynom
    Pn = legendre(l,cos(theta),'norm');
    Pn = Pn(abs(m)+1,:);
    % sign change required for odd positive m
    if m >= 0 && mod(m, 2) == 1
      Pn = -Pn;
    end
    
    % compute the normalizer
    if(isempty(shouldNormalize))
        % the default is yes, we should normalize
        shouldNormalize = true;
    end
    if(~shouldNormalize)
        normalizer = 1;
    else
        normalizer = (2*l + 1) / (4*pi);
        normalizer = normalizer * quotientOfFactorials(l-m, l+m);
        normalizer = sqrt(normalizer);
    end
    
    % compute the real part of the exp
    if(m > 0)
        E = cos(m * phi);
    elseif(m < 0)
        E = sin(m * phi);
    else
        E = cos(0 * phi);
    end
    
    % finish
    Ylm = normalizer * Pn .* E;
    Ylm = Ylm';
end

function q = quotientOfFactorials(numerator_arg, denominator_arg)
    % @return: (numerator_arg)! / (denominator_arg)!
%     numerator_v = 1:numerator_arg;
%     denominator_v = 1:denominator_arg;
    if numerator_arg > denominator_arg
        numerator_v = denominator_arg + 1:numerator_arg;
        denominator_v = 1;
    else
        numerator_v = 1;
        denominator_v = numerator_arg + 1:denominator_arg;
    end
    q = prod(numerator_v) / prod(denominator_v);
end
