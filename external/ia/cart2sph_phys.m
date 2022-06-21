% Converts Cartesian coordinates to physics-convention spherical
% coordinates, and vice versa.
% 
% angles = cart2sph_phys(xyz, rev)
% 
% xyz:     N x 3 matrix of Cartesian coordinates of N points.
% rev:     If true, reverses the direction of the conversion; i.e. xyz and
%          angles will be the spherical and Cartesian coordinates,
%          respectively (default: false).
% 
% angles:  N x 3 matrix of spherical coordinates of the points.
%
% See also:  reconCSAODF, reconCSAODF3Q, sampleODFs, CSAODF_CLI, EXAMPLE,
%            EXAMPLE_CLI.
 
% Code by Iman Aganj.

function angles = cart2sph_phys(xyz, rev)

if exist('rev', 'var') && rev  % sph2cart_phys
    angles = [sin(xyz(:,1)).*cos(xyz(:,2)), sin(xyz(:,1)).*sin(xyz(:,2)), cos(xyz(:,1))];
    if size(xyz,2)>2
        angles = bsxfun(@times, xyz(:,3), angles);
    end
else
    [TH, PHI, R] = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));
    angles = [pi/2-PHI, TH, R];
end
