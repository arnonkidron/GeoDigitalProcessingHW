
dir = "C:\Users\Arnon\Documents\GitHub\GeoDigitalProcessingHW\OFF models\";
mesh = MeshBasic(dir + "sphere_s1.off");
sphereCoords = getSphericalCoordinates(mesh.Vertices);

max_degree = 4;
A = ones(mesh.numV, (max_degree+1) ^2);
for l=0:max_degree
    for m=0:l
        theta = sphereCoords(:,1);
        phi = sphereCoords(:,2);
        harmonic = Ylm(l,m,theta,phi, false);
        
        index = m + l*(max_degree+1) + 1;
        A(:, index) = harmonic;
    end
    for m=l+1:max_degree
        index = m + l*(max_degree+1) + 1;
        A(:,index) = NaN;
    end
end

RenderSeveralFunctions(mesh, A, []);

function RenderSphericalHarmonic(mesh, sphereCoords, l, m)
    theta = sphereCoords(:,1);
    phi = sphereCoords(:,2);
    harmonic = Ylm(l,m,theta,phi);
%     rgb2 = [(real(harmonic) + 1) / 2, zeros(mesh.numV, 1), (imag(harmonic) + 1) / 2];
%     
%     rgb = mat2rgbCplx(harmonic, 1, 1);
%     rgb = reshape(rgb,[mesh.numV, 3]);
%     %Render(mesh, rgb);
    
    n = vecnorm(harmonic,2,2);
    Render(mesh, n);
    
end

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
