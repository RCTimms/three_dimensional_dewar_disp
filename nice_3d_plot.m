clc; clearvars; close all
cd /Users/rtimms/Documents/DPhil/man_dec_2020

set(0, 'DefaultAxesFontName', 'Helvetica');
set(0, 'DefaultTextFontName', 'Helvetica');

addpath(['/Users/rtimms/Documents/Source_recon_project/'...
    'MATLAB_concept_visualisation/npy-matlab-master/'])
addpath('/Users/rtimms/Downloads/osl/brewer_maps')


D_coreg=spm_eeg_load( '/Volumes/TASER/Notts_motor/03677_processed/03677341_Ellie_20170615_04/edff03677341_Ellie_20170615_04.mat')


cd(D_coreg.path)
% Load the gain mat
load(D_coreg.inv{1}.gainmat)

% Get the sensor info
sens=D_coreg.sensors('MEG');


% Get a FT vol object
vol=D_coreg.inv{1}.forward.vol;

% Channels kept in the lead field
D_coreg.inv{1}.inverse.Ic{1}

% % Show the sensors



clc;close all
sens.chanpos(find(strcmp(sens.chantype,'meggrad')),:)


[k,av] = convhull(ans,'simplify',true);
figure;plot3(ans(k,1),ans(k,2),ans(k,3))

p=trisurf(k,ans(:,1),ans(:,2),ans(:,3),'FaceColor','c');
mesh=[];
mesh.tri = p.Faces;
% For now, manually remove the obviously wrong faces. This could be
% automated by checking for euclidian distances between points, i.e. faces
% which have huge surface areas are likely to be errors
mesh.tri([251,252,254,255,77,253,188,417,418,315],:)=[];
mesh.pos = p.Vertices;

% Could use something like this, for example
% verts = get(p, 'Vertices');
% faces = get(p, 'Faces');
% a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
% b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));
% fprintf('\nThe surface area is %f\n\n', area);

close all;




figure
hold all;
% scatter3(ans(:,1),ans(:,2),ans(:,3),'ro','filled')
ft_plot_sens(sens,'facecolor','none', 'coilshape','point','style','k')

for i=350
data=D_coreg(3:273,i,1);

spm_mesh=[];
spm_mesh.faces=mesh.tri;
spm_mesh.vertices=mesh.pos;
smoothed_data=spm_mesh_smooth(spm_mesh,data,10)

% data=normalize(data,'range');


ft_plot_mesh(mesh,'vertexcolor',smoothed_data,'edgecolor','none','faceindex','False','facealpha',0.81)

colormap jet
drawnow
end


colormap(brewermap(128,'RdBu'))
spm_eeg_inv_checkdatareg_3Donly(D_coreg)

Surf = {mesh.tri,[ans(:,1), ans(:,2) ans(:,3)]};   
% IsoLine(Surf,data,10)

%%
% EXAMPLE 1: Isolines of f(x,y,z) = xy+20z on a triangular mesh -----------
figure(1), clf reset                       
set(1,'Name','Isolines of f(x,y,z) = xy+20z on a triangular mesh')
load trimesh3d                           % sample triangulation
F = x.*y+ 20*z;                          % function values 
trisurf(tri,x,y,z,F), shading interp     % plot surface, colored by F
Surf = {tri,[x y z]};                    % surface as cell array
[~,V]= IsoLine(Surf,F);                  % plot isolines        
camlight, axis equal                     % figure settings
zoom(1.3),view(-35,-27)
%
% Basically, that's all - but let's add some gimmicks ...
%
% Draw a colorbar, labelled with the iso-values
J = (.5+hsv(21))/1.5;                    % colormap, slightly whitened
J = J(floor((3:65)/3),:);                % triple entries
colormap(J);                             % apply J
c = colorbar('YTick',round(V));

function [H,V] = IsoLine(Surf,F,V,Col)
% IsoLine    plot isolines of a function on a surface
%    Unlike programs from the CONTOUR-family, IsoLine can handle 
%    arbitrary functions defined on the surface. Further, the 
%    surface can be given as a rectangular or a triangular mesh.
%
%    The call [H,V] = IsoLine(Surf,F,V,Col) is defined as follows:
%
%    Surf is a cell array containing a surface representation. It is
%      * either a triangular mesh {T,[x y z]}, where the rows of T  
%      are indices of vertices forming triangles, and the rows of
%      [x y z] are the coordinates of the vertices,
%      * or a rectangular mesh {X,Y,Z} with matrices of equal size.
%    F are the function values at the vertices, given as a vector 
%      the size of x or as a matrix the size of X, respectively.
%    V determines the values of F to be traced (default V=20).
%      As in CONTOUR, a singleton specifies the number of equidistant 
%      values spread over the range of F, while a vector specifies 
%      specific values. Use V = [v v] to draw a single line at level v.
%    Col is a character defining the color of isolines (default Col='k').
%
%    H is a vector containing the grahics handles of the isolines.
%    V is a vector containing the traced function values, what may be 
%      useful if V was given as a number.
%
% See also: IsoLineDemo, contour, contour3, delaunay, mesh
% -------------------------------------------------------------------------
% Author:  Ulrich Reif
% Date:    March 12, 2015
% Version: Tested with MATLAB R2014b and R2012b
% -------------------------------------------------------------------------
% Preprocess input --------------------------------------------------------
if length(Surf)==3                % convert mesh to triangulation
  P = [Surf{1}(:) Surf{2}(:) Surf{3}(:)];
  Surf{1}(end,:) = 1i;
  Surf{1}(:,end) = 1i;
  i = find(~imag(Surf{1}(:)));
  n = size(Surf{1},1);
  T = [i i+1 i+n; i+1 i+n+1 i+n];
else
  T = Surf{1};
  P = Surf{2};
end
f = F(T(:));
if nargin==2
  V = linspace(min(f),max(f),22);
  V = V(2:end-1);
elseif numel(V)==1
  V = linspace(min(f),max(f),V+2);
  V = V(2:end-1);
end
if nargin<4
  Col = 'k';
end
H = NaN + V(:);
q = [1:3 1:3];
% -------------------------------------------------------------------------
% Loop over iso-values ----------------------------------------------------
for k = 1:numel(V)
  R = {[],[]};
  G = F(T) - V(k);   
  C = 1./(1-G./G(:,[2 3 1]));
  f = unique(T(~isfinite(C))); % remove degeneracies by random perturbation
  F(f) = F(f).*(1+1e-12*rand(size(F(f)))) + 1e-12*rand(size(F(f)));
  G = F(T) - V(k);
  C = 1./(1-G./G(:,[2 3 1]));
  C(C<0|C>1) = -1;
  % process active triangles
  for i = 1:3
    f = any(C>=0,2) & C(:,i)<0;
    for j = i+1:i+2
      w = C(f,q([j j j]));
      R{j-i} = [R{j-i}; w.*P(T(f,q(j)),:)+(1-w).*P(T(f,q(j+1)),:)];
    end
  end
  % define isoline
  for i = 1:3
    X{i} = [R{1}(:,i) R{2}(:,i) nan+R{1}(:,i)]';
    X{i} = X{i}(:)';
  end
  % plot isoline
  if ~isempty(R{1})
    hold on
    H(k) = plot3(X{1},X{2},X{3},Col);
  end
end
% -------------------------------------------------------------------------
% Render with 'zbuffer' for best results ----------------------------------
set(gcf,'Renderer','zbuffer')
% -------------------------------------------------------------------------
end