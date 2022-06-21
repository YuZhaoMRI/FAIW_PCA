% Creates the matrix 'T' (and its pseudo-inverse 'pT') that projects the
% directions specified in 'angles' to a real and symmetric spherical
% harmonic basis of the order 'basisOrder'. The regularization 'lambda' in
% (Descoteaux et al, MRM'07) can optionally be specified, in which case the
% Laplacian 'L' (that can be made using makeL.m) can be optionally given.
% 
% [T, pT] = makePT(basisOrder, angles, lambda, L)
%
% See also:  compSH, cart2sph_phys, makeL, reconCSAODF, reconCSAODF3Q,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI. 

% Code by Iman Aganj.

function [T, pT] = makePT(basisOrder, angles, lambda, L)

T = compSH(basisOrder, angles(:,1)', angles(:,2)')';

if nargout>1
    if exist('lambda', 'var')
        if ~exist('L', 'var')
            L = makeL(basisOrder);
        end
        pT = (T'*T + lambda*(L'*L)) \ T';
    else
        pT = pinv(T);
    end
end
