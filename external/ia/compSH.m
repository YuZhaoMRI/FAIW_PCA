% Computes the spherical harmonic functions 'Y' of order 'basisOrder' and
% their partial derivatives 'dY' along the directions specified by row
% vectors 'theta' and 'phi'. 'precompF' (made by precompSH.m), 'cTheta'
% (cos(theta)), and 'sTheta' (sin(theta)) can be pre-computed and provided
% for speed-up.
% 
% [Y, dY] = compSH(basisOrder, theta, phi, precompF, cTheta, sTheta)
%
% See also:  precompSH, legendreSH, makePT, cart2sph_phys, reconCSAODF,
%            reconCSAODF3Q, sampleODFs, showODFs, CSAODF_CLI, EXAMPLE,
%            EXAMPLE_CLI.

% Code by Iman Aganj.

function [Y, dY] = compSH(basisOrder, theta, phi, precompF, cTheta, sTheta)

if nargin < 4
    precompF = precompSH(basisOrder);
end
if nargin < 5
    cTheta = cos(theta);
    sTheta = sin(theta);
end
oldVer = verLessThan('matlab', '9.1'); % Matlab releases earlier than R2016b.

nBasisFunctions = (basisOrder+1)*(basisOrder+2)/2;
nAngles = length(phi);
Y = zeros(nBasisFunctions, nAngles, 'like', theta);
for k = 0:2:basisOrder
    Pkm = legendreSH(k, cTheta, sTheta);
    for m = -k:k
        j = k*(k+1)/2+m+1;
        if m<0
            Y(j,:) = Pkm(-m+1,:).*cos(m*phi);
        elseif m==0
            Y(j,:) = Pkm(1,:);
        else
            Y(j,:) = Pkm(m+1,:).*sin(m*phi);
        end
        %j2 = k*(k-1)/2 + 2*abs(m) - (m<0) + 1;
    end
end
if ~isempty(precompF)
    if ~oldVer || nAngles == 1
        Y = Y .* precompF;
    else
        Y = bsxfun(@times, Y, precompF);
    end
end

if nargout > 1
    dY = zeros(nBasisFunctions, nAngles, 2, 'like', theta);
    for k = 0:2:basisOrder
        [Pkm, dPkm] = legendreSH(k, cTheta, sTheta);
        if isempty(dPkm)
            dY = [];
            return
        end
        for m = -k:k
            j = k*(k+1)/2+m+1;
            if m<0
                dY(j,:,1) = dPkm(-m+1,:).*cos(m*phi);
                dY(j,:,2) = -m*(Pkm(-m+1,:).*sin(m*phi));
            elseif m==0
                dY(j,:,1) = dPkm(1,:);
                dY(j,:,2) = 0;
            else
                dY(j,:,1) = dPkm(m+1,:).*sin(m*phi);
                dY(j,:,2) = m*(Pkm(m+1,:).*cos(m*phi));
            end
        end
    end
    if ~isempty(precompF)
        if oldVer
            dY = bsxfun(@times, dY, precompF);
        else
            dY = dY .* precompF;
        end
    end
end
