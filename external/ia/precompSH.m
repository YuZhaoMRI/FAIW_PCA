% Pre-computes some of the factors needed by compSH, for speed-up.
% 
% factors = precompSH(basisOrder)
%
% See also:  compSH, legendreSH, reconCSAODF, reconCSAODF3Q, CSAODF_CLI,
%            EXAMPLE, EXAMPLE_CLI.

% Code by Iman Aganj.

function factors = precompSH(basisOrder)

factors = zeros((basisOrder+1)*(basisOrder+2)/2, 1, 'like', basisOrder);
for k = 0:2:basisOrder
    for m = -k:k
        j = k*(k+1)/2+m+1;
        if m<0
            factors(j) = sqrt(((2*k+1)/(2*pi))*factorial(k+m)/factorial(k-m));
        elseif m==0
            factors(j) = sqrt((2*k+1)/(4*pi));
        else
            factors(j) = (-1)^m*sqrt(((2*k+1)/(2*pi))*factorial(k-m)/factorial(k+m));
        end
        %j2 = k*(k-1)/2 + 2*abs(m) - (m<0) + 1;
    end
end
