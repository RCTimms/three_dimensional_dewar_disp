function show_data_on_dewar(D_coreg,data_to_show,show_iso_lines,disp_anat,N_isolines,tpts_of_interest,col_map);
% Need to add documentation to inputs and assign defaults here:


% Get the sensor info. Should be capable of plotting both grad and ref data
% in the future, and dealing with both inverted and non-inverted D objects.
% This being said, the D objects will always have to be co-regged.

% Ensure data are of the right type
data_to_show=double(data_to_show);

% Deal with bad channels
good_meeg_chans=D_coreg.indchantype('MEEG','GOOD');
vol=D_coreg.inv{1}.forward.vol;
sens=D_coreg.inv{1}.forward.sensors;

% Prepare the sensor object with only the remaining good channels. This is
% particularly important when using SPM's inverse methods as any data
% dimensionality will **only be applied to the good channels!**
[~,sens]=ft_prepare_vol_sens(vol,sens,'channel',D_coreg.chanlabels(good_meeg_chans));


meg_sensors=sens.chanpos(find(strcmp(sens.chantype,'meggrad')),:);

% Approximate the dewar shape with the convhull function
[k,av] = convhull(meg_sensors,'simplify',true);
p=trisurf(k,meg_sensors(:,1),meg_sensors(:,2),meg_sensors(:,3),'FaceColor','c');
mesh=[];
mesh.tri = p.Faces;

% For now, manually remove the obviously wrong faces. This could be
% automated by checking for euclidian distances between points, i.e. faces
% which have huge surface areas are likely to be convhull "errors"
% mesh.tri([251,252,254,255,77,253,188,417,418,315],:)=[];
mesh.pos = p.Vertices;
close all; % close the temporary figure (annoying)


% We could use something like this, to work out the surface area of the
% faces, for example
% verts = get(p, 'Vertices');
% faces = get(p, 'Faces');
% a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
% b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));
% fprintf('\nThe surface area is %f\n\n', area);



%%% Make plot
figure
hold all;
ft_plot_sens(sens,'facecolor','none', 'coilshape','point','style','k')

for i=tpts_of_interest % Time point of interest
    
    % Just get the MEG data - should code this with something like
    % strcmp('MEGGRAD')... etc
    data=data_to_show;
    
    % Smooth the data with SPM. 5 iterations seems about appropriate. Can
    % try and remove this dependency in the future
    spm_mesh=[];
    spm_mesh.faces=mesh.tri;
    spm_mesh.vertices=mesh.pos;
    try
        smoothed_data=spm_mesh_smooth(spm_mesh,data,5);
    catch
        warning('Could not smooth data with spm_mesh_smooth. Data displayed will be the original');
        smoothed_data=data;
    end
    
    % Plot the mesh
    ft_plot_mesh(mesh,'vertexcolor',smoothed_data,'edgecolor','none','faceindex','False','facealpha',0.81)
    
    % Make the colourmap nice and sexy, if possible
    try
        colormap(brewermap(128,col_map))
    catch
        ft_hastoolbox('BREWERMAP',0,0)
        colormap(jet)
    end
    drawnow
end
axis on
pause(2)
% The following line will plot the subject anatomy within the dewar. We can
% extract the code from this function and make it pure FieldTrip with
% relative ease in the future
if disp_anat==1;hold all;spm_eeg_inv_checkdatareg_3Donly(D_coreg);end

% Show magnetic field isolines?
if show_iso_lines==1;
    Surf = {mesh.tri,[meg_sensors(:,1), meg_sensors(:,2) meg_sensors(:,3)]};
    IsoLine(Surf,data,N_isolines);
end
end


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