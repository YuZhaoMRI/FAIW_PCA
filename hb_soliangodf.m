function  [odfssi,Nt] = hb_soliangodf(odfs,verts,ndirs,omega,shord,nico,sumint,verif)
%HB_SOLIANGODF computes sum of samples (or numerical integral) of odfs
% within constant solid angles subtended along a given set of directions. 
% The odfs are first represented in a spherical harmonic basis and then 
% uniformly resampled within the solid angles. The sampling points are 
% determined using the vertices of a subdivided icosahedron. 
%
% Inputs:
%   odfs: MxN matrix of odf samples; one odf per column.
%   vertices: 3xM matrix of unit norm odf vertices.
%   ndirs: 3xD matrix of unit norm neighborhood directions.
%   omega: (optional) solid angle. default: 4*pi/size(ndirs,2)
%   shord: (optional) spherical harmonic basis order
%   nico: (optional) number of icosahedron divisions. default: 5
%   sumint: (optional) 'sum' or 'integral'. default: 'sum' 
%   verif: (optional) verify template & rotations. default: false
%
% Outputs:
%   odfssi: sum/integral of odfs within solid angles. 

% Hamid Behjat
% May 2020.

D = size(ndirs,2);
N = size(odfs,2);
if ~exist('shord','var')||isempty(shord)
    shord = 8;
end
if ~exist('omega','var')||isempty(omega)
    omega = 4*pi/D;
end
if ~exist('nico','var')||isempty(nico)
    nico = 5;
end
if ~exist('sumint','var')||isempty(sumint)
    sumint = 'sum';
end
if ~exist('verif','var')||isempty(verif)
    verif = false;
end

% odfs in sh domain
angles = cart2sph_phys(verts');
[~,pT] = makePT(shord,angles);
sh = pT*odfs;

% build sampling template
costheta = 1-omega/(2*pi);
d = icosphere(nico);
verts_icosph = d.Vertices';
zdir = [0 0 1]';
d1 = zdir'*verts_icosph>costheta;
templ = verts_icosph(:,d1);
if verif,hax = plotsamps('templ',templ);end
centind = find(all(eq(zdir,templ)));
Nt = size(templ,2);

odfssi = zeros(D,N);
for iD=1:D
    d = ndirs(:,iD);
    % rotate template
    R = rotmatrix(zdir,d);
    d1 = R*templ; 
    if verif,plotsamps('ndir',d1,d,iD,hax,centind);end
    % sample odfs
    d2 = cart2sph_phys(d1');
    T = makePT(shord,d2);
    % sum
    odfssi(iD,:) = sum(T*sh);
end
if strcmp(sumint,'integral')
    % dS: omega/Nt
    odfssi = odfssi*(omega/Nt);
end
end

function R=rotmatrix(v1,v2)
% Rotation matrix tmapping v1 to v2: v2=R*v1
d = v1'*v2;
c = cross(v1,v2);
nb = norm(v1);
nd = norm(v2);
nc = norm(c);
if ~all(c==0)
    z = [0 -c(3) c(2); c(3) 0 -c(1); -c(2) c(1) 0];
    R = (eye(3)+z+z^2*(1-d)/nc^2)/nb^2;
else % collinear
    R = sign(d)*(nd/nb) ;
end
end

function ha=plotsamps(opts,varargin)
colorblind = true;
n = size(varargin{1},2);
d = hsv(n);
rng(1);
clmap = d(randperm(n),:);
switch opts
    case 'templ'
        d = findobj('type','figure');
        if isempty(d)
            figure(1);
        else
            figure(max(arrayfun(@(x) x.Number,d))+1);
        end
        [x,y,z] = sphere(50);
        r = .99;
        hS = surf(r*x,r*y,r*z,'FaceAlpha',1);
        set(hS,'FaceColor',.9*ones(1,3),'EdgeColor','none');
        hold on
        p = varargin{1};
        scatter3(p(1,:)',p(2,:)',p(3,:)',15,'k');
        text(p(1,end),p(2,end),p(3,end),' template','FontSize',13);
        set(gcf,'Position',get(0,'Screensize'));
        view(3)
        axis tight
        axis equal
        ha=gca;
    case 'ndir'
        p = varargin{1};
        l = 1.05*varargin{2};
        cl = clmap(varargin{3},:);
        ha = varargin{4};
        plot3(ha,[0,l(1)],[0,l(2)],[0,l(3)],'Color',cl,'LineWidth',3)
        scatter3(ha,p(1,:)',p(2,:)',p(3,:)',7,cl,'filled');
        % The color of the sticks & samples should match.
        % Can also be verified by checking if tags match:
        if colorblind
            t = num2str(varargin{3});
            s = 'HorizontalAlignment';
            ci = varargin{5};
            text(ha,l(1),l(2),l(3),t,'FontSize',10,s,'left');
            text(ha,p(1,ci),p(2,ci),p(3,ci),t,'FontSize',10,s,'right');
            %view(l')
        end
end
drawnow
figure(gcf)
pause(.01)
%pause
end