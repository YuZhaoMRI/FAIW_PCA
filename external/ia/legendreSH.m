% Similarly to Matlab's legendre function, legendreSH samples the
% associated Legendre polynomials 'Pkm' of degree 'k', however faster for
% values of k<=6. The directions must be specified by row vectors 'c'
% (cos(theta)) and 's' (sin(theta)). The derivatives 'dPkm' with respect to
% theta are also computed.
% 
% [Pkm, dPkm] = legendreSH(k, c, s)
%
% See also:  legendre, compSH, cart2sph_phys, reconCSAODF, reconCSAODF3Q,
%            CSAODF_CLI, EXAMPLE, EXAMPLE_CLI. 

% Code by Iman Aganj.

function [Pkm, dPkm] = legendreSH(k, c, s)

switch k
    case 0
        Pkm = ones(size(c), 'like', c);
    case 2
        Pkm = [1.5*c.^2-0.5; -3*c.*abs(s); 3*s.^2];
    case 4
        c2 = c.^2;
        s2 = 1-c2;
        cs = c.*abs(s);
        Pkm = [(4.375*c2-3.75).*c2+0.375; (7.5-17.5*c2).*cs; (52.5*c2-7.5).*s2; -105*cs.*s2; 105*s2.^2];
    case 6
        c2 = c.^2;
        c4 = c2.^2;
        c6 = c2.*c4;
        cs = c.*abs(s);
        Pkm = [14.4375*c6 - 19.6875*c4 + 6.5625*c2 - 0.3125; cs.*(-86.625*c4 + 78.75*c2 - 13.125); -433.125*c6 + 669.375*c4 - 249.375*c2 + 13.125; cs.*(1732.5*c4 - 2205*c2 + 472.5); 5197.5*c6 - 10867.5*c4 + 6142.5*c2 - 472.5; cs.*(-10395*c4 + 20790*c2 - 10395); -10395*c6 + 31185*c4 - 31185*c2 + 10395];
    otherwise
        Pkm = legendre(k, c);
end

if nargout > 1
    switch k
        case 0
            dPkm = zeros(size(c), 'like', c);
        case 2
            c2 = c.^2;
            s2 = 1-c2;
            cs = c.*s;
            dPkm = [-3*cs; -3*(c2-s2).*sign(s); 6*cs];
        case 4
            cs = c.*s;
            sgnS = sign(s);
            dPkm = [(7.5-17.5*c2).*cs; ((7.5-17.5*c2).*(c2-s2)+35*c2.*s2).*sgnS; (90-210*s2).*cs; 105*s2.*(1-4*c2).*sgnS; 210*cs];
        case 16
            cs = c.*s;
            sgnS = sign(s);
            dPkm = [(-86.625*c4 + 78.75*c2 - 13.125).*cs; (-519.75*c6 + 748.125*c4 - 262.5*c2 + 13.125).*sgnS; (2598.75*c4 + - 2677.5*c2 + 498.75).*cs; (10395*c6 - 17482.5*c4 +  7560*c2 - 472.5).*sgnS; (-31185*c4 + 43470*c2 -12285).*cs; (-62370*c6 + 135135*c4 - 83160*c2 + 10395).*sgnS; (62370*c4 - 124740*c2 + 62370).*cs];
        otherwise
            dPkm = bsxfun(@rdivide, k*bsxfun(@times, c, Pkm) - [bsxfun(@times, k+(0:k-1)', legendre(k-1, c)); zeros(1,length(c))], abs(s));
    end
end
