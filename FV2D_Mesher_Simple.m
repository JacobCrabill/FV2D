% -------------------------------------------------------------------------
%
%                          FV2D Mesh Generator
%   Creates geometry and data structures necessary for FV2D Euler solver
%                       
%                     AA 214C, Stanford University
%                       Jacob Crabill  4/21/2014
% -------------------------------------------------------------------------
clear all; close all;

%% ------------------------------------------------------------------------
% Mesh Parameters
% -------------------------------------------------------------------------
plotmesh = true;

% Mesh File Name
meshFile = 'mesh/testBox';

% Wall-Plus-Bump (true) or flow past wedge (false)
bumpMesh = false;

% Supersonic or no?
supersonic = false;

% Viscous or inviscid?
viscous = true;

% Top/Bottom surface is slip wall (true) or inlet/outlet (false)?
topWall = true;
btmWall = true;

% Bump, starting at (0,0)
L = 1;

% Overall dimensions
xmin = 0;
ymin = 0;
xmax = L;
ymax = L;

% Number of CELLS in mesh (NOT vertices)
nx = 5;
ny = 5;

%% ------------------------------------------------------------------------
% Mesh Generation
% -------------------------------------------------------------------------
% Begin with Cartesian grid, [xi,eta] in range [0,1]
% Scale accordingly to create stretched grid with bump
disp('  Beginning mesh generation...');

% ------------------------------------------
% Create the Mesh
% ------------------------------------------
% Create Cartesian domain with proper # of cells
dxi = 1/nx;
deta = 1/ny;
xi = 0:dxi:1;
eta = 0:deta:1;

expScale = 1.1;
eta = (eta-ymin).*(1/(1-exp(expScale))).*(1-exp(expScale*eta))+ymin;
[XI,ETA] = meshgrid(xi,eta);

X = XI*(xmax-xmin) + xmin;
Y = ETA*(ymax-ymin) + ymin;

% Show the mesh vertices
if plotmesh
    plot(X,Y,'k*');
end

%% ------------------------------------------------------------------------
% Setup Mesh Data Storage
% -------------------------------------------------------------------------
xv = zeros((nx+1)*(ny+1),2);
c2v = zeros(nx*ny,4);
c2nf = ones(nx*ny,1)*4;

% Setup Vertices
nv = 0;
for j=1:(ny+1)
    for i=1:(nx+1)
        nv = nv + 1;
        xv(nv,:) = [X(j,i),Y(j,i)];
    end
end

% Setup Cells
nc = 0;
for i=1:nx
    for j=1:ny
        nc = nc + 1;
        c2v(nc,1) = (j-1)*(nx+1)+i;
        c2v(nc,2) = (j-1)*(nx+1)+i+1;
        c2v(nc,3) = j*(nx+1)+i+1;
        c2v(nc,4) = j*(nx+1)+i;
    end
end
n_cells = size(c2v,1);

% ------------------------------------------
% Setup Boundaries, BC's
% ------------------------------------------
% 1 = subsonic inlet
% 2 = subsonic outlet
% 3 = supersonic inlet
% 4 = supersonic outlet
% 5 = slip wall
% 6 = isothermal no-slip wall
% 7 = adiabatic no-slip wall
n_bnds = 5;
bounds = cell(5);

bounds{1}.name = 'inlet';
bounds{2}.name = 'outlet';
bounds{3}.name = 'top';
bounds{4}.name = 'wall';
bounds{5}.name = 'symmetry'; % pre-leading edge btm srfc

if supersonic
    bounds{1}.tag = 3;  % Inlet
    bounds{2}.tag = 4;  % Outlet
    if topWall
        bounds{3}.tag = 5;  % Top / 'Freestream'
    else
        bounds{3}.tag = 3;  % Top / 'Freestream'
    end
else
    bounds{1}.tag = 1;  % Inlet
    bounds{2}.tag = 2;  % Outlet
    if topWall
        bounds{3}.tag = 5;  % Top / 'Freestream'
    else
        bounds{3}.tag = 1;  % Top / 'Freestream'
    end
end

% Bottom / Wall

if btmWall
    if viscous
        %btmBC = 6;  % Isothermal Wall
        btmBC = 7;  % Adiabatic Wall
    else
        btmBC = 5;  % Slip Wall
    end
    btmBC2 = 5;
else
    if supersonic
        btmBC = 4;
        btmBC2 = 4;
    else
        btmBC = 2;
        btmBC2 = 2;
    end
end

bounds{4}.tag = btmBC;  % adiabatic wall
bounds{5}.tag = btmBC2; % slip wall / symmetry plane

bounds{1}.pts = 1:(nx+1):((nx+1)*ny+1);
bounds{2}.pts = (nx+1):(nx+1):((nx+1)*(ny+1));
bounds{3}.pts = (1:(nx+1))+ny*(nx+1);
bounds{4}.pts = 1:nx+1;

%% ------------------------------------------------------------------------
% Mesh Pre-Processing
% -------------------------------------------------------------------------
% To be done on any mesh that comes in (not just generated in-place)
% ------------------------------------------
% Setup Edges
% ------------------------------------------
max_nf = max(c2nf);
e2v = zeros(nx*ny*max_nf,2);
e2c1 = -1*ones(nx*ny*max_nf,1);
c2e1 = -1*ones(nx*ny*max_nf,1);
c2e = -1*ones(nx*ny*max_nf,2);

ne = 0;
for ic = 1:n_cells
    for k=1:c2nf(ic)
        if k == c2nf(ic)
            kp = 1;
        else
            kp = k + 1;
        end
        ne = ne + 1;
        e2v(ne,:) = c2v(ic,[k,kp]);
        e2c1(ne) = ic;
        c2e1(ic,k) = ne;
    end
end

e2v1 = sort(e2v,2);
% ia = indicies of the first unique row in e2v1
% iE = indicies of e2v which map to each row of e2v1
[e2v,~,iE] = unique(e2v1,'rows');

% ------------------------------------------
% Setup Internal & Boundary Edge Lists
% ------------------------------------------
e2c = zeros(size(e2v));
bnde = zeros(length(e2v),1);  % boundary interfaces
inte = zeros(length(e2v),1);  % internal interfaces
bc = zeros(length(e2v),1);    % boundary conditions
n_bnd_edges = 0; n_int_edges = 0; n_edges = length(e2v);
for i=1:size(iE,1)
    if iE(i) ~= -1
        ie = find(iE==iE(i));
        if length(ie) > 2
            error(['More than 2 cells for edge ',i]);
        elseif length(ie) == 2
            % internal edge which has not yet been added
            %e2c(iE(i),:) = [e2c1(ie(1)),e2c1(ie(2))];
            n_int_edges = n_int_edges + 1;
            inte(n_int_edges) = iE(i);
        elseif length(ie) == 1
            % boundary edge
            for k = 1:n_bnds
                ib1 = find(bounds{k}.pts == e2v(iE(i),1));
                ib2 = find(bounds{k}.pts == e2v(iE(i),2));
                if numel(ib1) > 0 && numel(ib2) > 0
                    n_bnd_edges = n_bnd_edges + 1;
                    bnde(n_bnd_edges) = iE(i);
                    bc(n_bnd_edges) = bounds{k}.tag;
                    break;
                end
            end
        end
        iE(iE==iE(i)) = -1;
    end
end
% adding in ghost cells
n_ghost_cells = n_bnd_edges;
n_cells_total = n_cells+n_ghost_cells;

% ------------------------------------------
% Setup Edge Normals
% ------------------------------------------
unorm = zeros(n_edges,2);
dA = zeros(n_edges,1);
for ie=1:n_edges
    dx = xv(e2v(ie,2),1) - xv(e2v(ie,1),1);
    dy = xv(e2v(ie,2),2) - xv(e2v(ie,1),2);
    theta = atan2(dy,dx);
    S = sin(theta); C = cos(theta);
    dA(ie) = sqrt(dx^2+dy^2);
    unorm(ie,:) = [S,-C]; % n_hat to 'right' of edge looking from 1-2
end

% ------------------------------------------
% Setup Cell-To-Edge
% ------------------------------------------
% Also Match cell OUTWARD normal direction to edge normal direction (+/- 1)
% Since edge normal points to 'right' of edge, outward normal corresponds
% to c2v(ic,i) < c2v(ic,i+1)  [CCW ordering]
c2n = zeros(n_cells_total,max_nf);
for ic = 1:n_cells
    for j = 1:max_nf
        if j>c2nf(ic)
            %c2e(ic,j) = n_edges+1; % this ain't gonna work
            error('triangles not yet working, sorry!');
        else
            jp1 = j+1;
            if j==c2nf(ic), jp1=1; end;

            if c2v(ic,j) < c2v(ic,jp1)
                c2n(ic,j) = 1;
                edge = [c2v(ic,j),c2v(ic,jp1)];
            else
                c2n(ic,j) = -1;
                edge = [c2v(ic,jp1),c2v(ic,j)];
            end

            ie1 = find(e2v(:,1)==edge(1));
            ie2 = find(e2v(ie1,2)==edge(2));
            if length(ie2) > 1
                error('Edge appears more than twice in e2v!');
            elseif numel(ie2) == 0
                error('Cell edge not found in e2v!');
            end
            ie = ie1(ie2);
            if ie < 0
                error('invalid edge number!');
            end
            c2e(ic,j) = ie;

            if c2n(ic,j) > 0
                e2c(ie,1) = ic;
            else
                e2c(ie,2) = ic;
            end
        end
    end
end

% ------------------------------------------
% Cleanup - Remove Padding from all arrays
% ------------------------------------------
% c2n, unorm sizes set properly to begin with
if length(c2e) > n_cells
    c2e(n_cells+1:length(c2e),:) = [];
end
if length(c2v) > n_cells
    c2v(n_cells+1:length(c2v),:) = [];
end
if length(e2c) > n_edges
    e2c(n_edges+1:length(e2c),:) = [];
end
if length(e2v) > n_edges
    e2v(n_edges+1:length(e2v),:) = [];
end
if length(bnde) > n_bnd_edges
    bnde(n_bnd_edges+1:length(bnde),:) = [];
end
if length(inte) > n_int_edges
    inte(n_int_edges+1:length(inte),:) = [];
end

% Throw in ghost cells, as well, to c2e (remaining 'ghost edges' to be
% added later)
c2e = [c2e; zeros(n_ghost_cells,4)];
c2e((n_cells+1):n_cells_total,1) = bnde;

% ------------------------------------------
% Finalize Boundary Edges
% ------------------------------------------
% Set ghost-cell index, and set interior cell to be on "left" of edge. Also
% adjust c2n accordingly. If g.c. already on right, just assign it an ID.
for ie = 1:n_bnd_edges
    if e2c(bnde(ie),1) == 0
        e2c(bnde(ie),:) = fliplr(e2c(bnde(ie),:)); % set g.c. to right side
        e2c(bnde(ie),2) = n_cells+ie;                  % assign g.c. index
        unorm(bnde(ie),:) = - unorm(bnde(ie),:);       % flip edge normal
        ic = e2c(bnde(ie),1);
        lfid = find(c2e(ic,:)==bnde(ie));    % find local face id within cell
        if isempty(lfid) 
            error('edge touches cell but cell does not have that edge!')
        end
        c2n(ic,lfid) = -c2n(ic,lfid);  % flip cell normal
    else
        e2c(bnde(ie),2) = n_cells+ie;  % assign ghost cell index
    end
	c2n(e2c(bnde(ie),2),1) = -1;    % reverse g.c. normal (since pointing into cell)
end

% ------------------------------------------
% Setup Ghost Cells
% ------------------------------------------
% Ghost cells always on 'right' of boundary edges
e2c(bnde,2) = n_cells + (1:n_bnd_edges)';
gc2v = zeros(n_ghost_cells,4);
gc2v(:,1:2) = e2v(bnde,:);

% ------------------------------------------
% Setup Cell-To-Cell
% ------------------------------------------
c2c = zeros(n_cells_total,max_nf);

ic = (1:n_cells)';
for j = 1:max_nf
    % Get the neighboring cell (index of each row of e2c != each ic)
    ICN = e2c(c2e(ic,j),:);
    ICN1 = find((ICN(:,1)-ic)~=0);
    ICN2 = find((ICN(:,2)-ic)~=0);
    ic2 = zeros(size(ic));
    ic2(ICN1) = ICN(ICN1,1);
    ic2(ICN2) = ICN(ICN2,2);

    ic1 = ic(ic2~=0);  % 'main' cell  (cell "i")
    ic2 = ic2(ic2~=0); % neighboring cell  (cell "j")
    
    c2c(ic1,j) = ic2;
end

% Take care of ghost cells - find interior cell neighbor.  This will be the
% first cell listed in 'c2c'
for ic=1:n_ghost_cells
    ic1 = e2c(bnde(ic),:);
    ic1 = ic1(ic1~=(ic+n_cells));
    c2c(ic+n_cells,1) = ic1;
end

% Take care of neighboring ghost cells (each g.c. has 1 or 2 neighbors).
% These will be the 2nd and 3rd entried in 'c2c'
c2nf = [c2nf; ones(n_ghost_cells,1)];
gc2v = gc2v(:,1:2);
for ic=1:n_ghost_cells
    [ic1,~] = find(gc2v==gc2v(ic,1));
    [ic2,~] = find(gc2v==gc2v(ic,2));
    if length(ic1) == 2
        ic1 = ic1(ic1~=ic);
        if c2c(ic+n_cells,1) ~= c2c(ic1+n_cells,1)
            c2nf(n_cells+ic) = c2nf(n_cells+ic) + 1;
            c2c(n_cells+ic,c2nf(n_cells+ic)) = ic1+n_cells;
        end
    end
    if length(ic2) == 2
        ic2 = ic2(ic2~=ic);
        if c2c(ic+n_cells,1) ~= c2c(ic2+n_cells,1)
            c2nf(n_cells+ic) = c2nf(n_cells+ic) + 1;
            c2c(n_cells+ic,c2nf(n_cells+ic)) = ic2+n_cells;
        end
    end
end

% Setup remaining ghost cell vertices
% Extrapolate along edge normal by distance equal to interior cell's width
%   2               1
%   X ============= X  <-- existing boundary edge & points
%   |               | 
%   |               | 
%   |       o       | 
%   |               | 
%   |               | 
%   x ------------- x  <-- ghost points created through extrapolation
%   3               4
gc2v = [gc2v,zeros(size(gc2v))];
gcxv = zeros(n_ghost_cells*2,2);
nv = 0;
for ic = 1:n_ghost_cells
    icn = c2c(n_cells+ic,1); % neighboring interior cell
    % find out which face of interior cell matches ghost cell
    i_face = find(c2c(icn,:) == n_cells+ic);
    if i_face == c2nf(icn)
        jfp = 1;
    else
        jfp = i_face + 1;
    end
    if i_face == 1
        jfm = c2nf(icn);
    else
        jfm = i_face - 1;
    end
    
    % see if existing nodes (1-2) are CW or CCW as seen from g.c. center
    gxv1 = xv(e2v(bnde(ic),1),:);
    gxv2 = xv(e2v(bnde(ic),2),:);
    
    xyc = (gxv1 + gxv1)/2 + unorm(bnde(ic),:);
    dx1 = [gxv1 - xyc, 0];
    dx2 = [gxv1 - xyc, 0];
    if cross(dx1,dx2) > 0   % already in CCW order
        gc2v(ic,3) = nv+1;
        gcxv(nv+1,:) = gxv2 + unorm(bnde(ic),:)*dA(c2e(icn,jfm));
        gc2v(ic,4) = nv+2;        
        gcxv(nv+2,:) = gxv1 + unorm(bnde(ic),:)*dA(c2e(icn,jfp));
    else                    % not in CCW order - flip
        gc2v(ic,1:2) = fliplr(gc2v(ic,1:2));
        gc2v(ic,3) = nv+1;
        gcxv(nv+1,:) = gxv1 + unorm(bnde(ic),:)*dA(c2e(icn,jfm));
        gc2v(ic,4) = nv+2;        
        gcxv(nv+2,:) = gxv2 + unorm(bnde(ic),:)*dA(c2e(icn,jfp));
    end 
    nv = nv + 2;    
end

% Merge ghost-cell vertices which should be shared between cells
GXV = inf*ones(size(gcxv));
gvtag = zeros(n_ghost_cells,2);
nv = 0;
for ic = 1:n_ghost_cells
    gv1 = gc2v(ic,1);
    gv2 = gc2v(ic,2);
    icn1 = find(gc2v(:,2)==gv1);
    icn2 = find(gc2v(:,1)==gv2);
    
    % for point 1, need to merge point 4 with point 3 of cell icn1
    % (If not done so already)
    % Also, is cell icn1 actually a neighbor, or is it around a corner?
    if ~isempty(icn1) && ~isempty(find(c2c(ic+n_cells,:)==icn1+n_cells, 1))
        if gvtag(icn1,1)<1 && gvtag(ic,2)<1
            gv4 = gc2v(ic,4);
            gv3 = gc2v(icn1,3);
            GXV(nv+1,:) = (gcxv(gv3,:)+gcxv(gv4,:))/2;
            gc2v(ic,4) = nv+1;
            gc2v(icn1,3) = nv+1;
            nv = nv + 1;
            % Mark the nodes as merged so we don't repeat the process
            gvtag(icn1,1) = 1;  % 1 / 2 corresponds to node 3 / 4 of cell
            gvtag(ic,2) = 1;            
        end
    else
        gv4 = gc2v(ic,4);
        GXV(nv+1,:) = gcxv(gv4,:);
        gc2v(ic,4) = nv+1;
        gvtag(ic,2) = 1;
        nv = nv + 1;
    end
    
    % for point 2, need to merge point 3 with point 4 of cell icn2
    % Also, is cell icn2 actually a neighbor, or is it around a corner?
    if ~isempty(icn2) && ~isempty(find(c2c(ic+n_cells,:)==icn2+n_cells, 1))
        if gvtag(icn2,2)<1 && gvtag(ic,1)<1
            gv4 = gc2v(icn2,4);
            gv3 = gc2v(ic,3);
            GXV(nv+1,:) = (gcxv(gv3,:)+gcxv(gv4,:))/2;
            gc2v(ic,3) = nv+1;
            gc2v(icn2,4) = nv+1;
            nv = nv + 1;
            % Mark the nodes as merged so we don't repeat the process
            gvtag(ic,1) = 1;
            gvtag(icn1,2) = 1; 
        end
    else
        % Keep the vertex as-is in the current ghost cell
        gv3 = gc2v(ic,3);
        GXV(nv+1,:) = gcxv(gv3,:);
        gc2v(ic,3) = nv+1;
        gvtag(ic,1) = 1;
        nv = nv + 1;
    end
end
% Convert ghost vertex indices to full-vertex-list indices
n_verts = length(xv);
gc2v(:,3:4) = gc2v(:,3:4) + n_verts;
% Trim down GXV to only final (merged) vertices
n_ghost_verts = nv;
n_verts_total = n_verts + n_ghost_verts;
GXV(nv+1:length(GXV),:) = [];

% ------------------------------------------
% Get list of ghost-cell edges & related info
% ------------------------------------------
ge2v = zeros(n_ghost_cells*3,2);
ge2c1 = zeros(n_ghost_cells*3,2);
gc2e1 = zeros(n_ghost_cells,3); % skipping bnde edges (already have 'em)

ne = 0;
for ic = 1:n_ghost_cells
    ne = ne + 1;
    ge2v(ne,:) = gc2v(ic,2:3);
    gc2e1(ic,1) = ne;
    
    ne = ne + 1;
    ge2v(ne,:) = gc2v(ic,3:4);
    gc2e1(ic,2) = ne;
    
    ne = ne + 1;
    ge2v(ne,:) = gc2v(ic,[4,1]);
    gc2e1(ic,3) = ne;    
end

% Remove padding from end
ge2v((ne+1):(n_ghost_cells*3),:) = [];

ge2v1 = sort(ge2v,2);
% ia = indicies of the first occurance of each unique row in ge2v1
% iE = indicies of ge2v which map to each row of ge2v1
[ge2v,ia,iE] = unique(ge2v1,'rows');

ic=(1:n_ghost_cells)';
gc2e = gc2e1;
for j=1:3
    gc2e(ic,j) = iE(gc2e1(ic,j));
end
gc2e = gc2e + n_edges; % convert to total-mesh indices
gc2e = [bnde,gc2e];    % combine known bnde edges + new ghost edges

% ------------------------------------------
% Merge ghost-cell vertices & edges with the rest
% ------------------------------------------
n_edges = length(e2v);
n_ghost_edges = length(ge2v);
n_edges_total = n_edges + n_ghost_edges;

xv = [xv;GXV];
c2v = [c2v;gc2v];
% ghost cells were actually part of c2e already, sort of... remove them!
c2e(n_cells+1:n_cells_total,:) = [];
c2e = [c2e;gc2e];
e2v = [e2v;ge2v];
c2nf(n_cells+1:n_cells_total) = 4;

% ------------------------------------------
% Vertices-to-cells (incl. ghost cells)
% ------------------------------------------
max_c_per_v = 20;
nv = length(xv);
v2c = zeros(nv,max_c_per_v);
v2nc = zeros(nv,1);
for i=1:nv
    [r,~] = find(c2v==i);
    %[rg,~] = find(gc2v(:,1:2)==i);
    %nr = length(r);
    %nrg = length(rg);
    
    v2c(i,1:length(r)) = r;
    %if nrg > 0
    %    rg = rg+n_cells;
    %    v2c(i,(nr+1):(nr+nrg)) = rg;
    %end
    %v2nc(i) = nr + nrg;
    v2nc(i) = length(r);
end
max_nc = max(v2nc);
v2c(:,max_nc+1:max_c_per_v) = [];

% ------------------------------------------
% Vertices-to-# of edges (incl. ghost edges?)
% ------------------------------------------
v2ne = zeros(n_verts_total,1);
for iv=1:n_verts_total
    v2ne(iv) = length(find(e2v==iv));
end

% ------------------------------------------
% Vertices to Edges, Vertices to Vertices
% ------------------------------------------
v2e = zeros(n_verts_total,max(v2ne));
v2v = zeros(n_verts_total,max(v2ne));
for iv=1:n_verts_total
    [ie1,~] = find(e2v==iv);
    v2e(iv,1:length(ie1)) = ie1;
    iv1 = e2v(ie1,:);
    iv1 = iv1(iv1~=iv);
    v2v(iv,1:length(ie1)) = iv1;
end

% ------------------------------------------
% Setup Edge-To-Local-Face
% ------------------------------------------
% Local face id (between 1:max_nf) for each of e2c
e2f = zeros(n_edges,2);
% e2f1 = zeros(n_edges,2);

ie = (1:n_edges)';
ic1 = e2c(ie,1);
ic2 = e2c(ie,2);

IE = repmat(ie,[1,max_nf]);

[r1,c1] = find(c2e(ic1,:)-IE == 0);
R1C1 = sortrows([r1,c1]);
[r2,c2] = find(c2e(ic2,:)-IE == 0);
R2C2 = sortrows([r2,c2]);

e2f(:,1) = R1C1(:,2);
e2f(:,2) = R2C2(:,2);

% % try 2
% for ic=1:n_cells_total
%     for j=1:c2nf(ic)
%         if c2e(ic,j)<=n_edges
%             if e2c(c2e(ic,j),1)==ic
%                 k = 1;
%             elseif e2c(c2e(ic,j),2)==ic
%                 k = 2;
%             else
%                 error('!!!');
%             end
%             e2f(c2e(ic,j),k) = j;
%         end
%     end
% end

% ------------------------------------------
% Separate each B.C. for vectorization
% ------------------------------------------
ie_sub_in = bnde(bc==1);
ie_sub_out = bnde(bc==2);
ie_sup_in = bnde(bc==3);
ie_sup_out = bnde(bc==4);
ie_slip_wall = bnde(bc==5);
ie_isothermal_wall = bnde(bc==6);
ie_adiabatic_wall = bnde(bc==7);

n_sub_in = length(ie_sub_in);
n_sub_out = length(ie_sub_out);
n_sup_in = length(ie_sup_in);
n_sup_out = length(ie_sup_out);
n_slip_wall = length(ie_slip_wall);
n_isothermal_wall = length(ie_isothermal_wall);
n_adiabatic_wall = length(ie_adiabatic_wall);

% ------------------------------------------
% Find Cell-Centers, Edge-Centers
% ------------------------------------------
% Ghost Cells
%gXYV = reshape([xv(gc2v(:,1:2),:);GXV(gc2v(:,3:4),:)],n_ghost_cells,4,2);
%gXYC = squeeze(mean(gXYV,2));
% center of each cell (including ghost cells)
XYV = reshape(xv(c2v,:),n_cells_total,max_nf,2);
XYC = squeeze(mean(XYV,2));
%XYC = [XYC;gXYC];
% center of each edge
XYV = reshape(xv(e2v,:),n_edges_total,2,2);
XYE = squeeze(mean(XYV,2));

%dXE = zeros(n_cells,max_nf,2);
dXE = zeros(2,max_nf,n_cells_total);
ic = (1:n_cells)';
for j = 1:max_nf
    dXE(:,j,:) = (XYE(c2e(:,j),:) - XYC(:,:))';
end

% ------------------------------------------
% Calculate Cell Areas
% ------------------------------------------
% Note: ghost cells a little different; nodes not necessarily ordered CCW
Area = zeros(n_cells_total,1);
for ic = 1:n_cells_total
    for j = 1:max_nf
        if j==1
            jm1 = c2nf(ic);
        else
            jm1 = j-1;
        end        
        if j>=c2nf(ic)
            jp1 = mod(j,c2nf(ic))+1;
        else
            jp1 = j+1;
        end
        Area(ic) = Area(ic) + xv(c2v(ic,j),1)*(xv(c2v(ic,jp1),2)-xv(c2v(ic,jm1),2));
    end
end
Area = abs(Area); % required thanks to ghost cell non-standard-ness

% ------------------------------------------
% Calculate dual-control-volume face areas
% ------------------------------------------
% Lengths of lines connecting cell centers across each edge
% ('real' edges only)
ic1 = e2c(:,1);
ic2 = e2c(:,2);
dA2 = sqrt(sum((XYC(ic2,:)-XYC(ic1,:)).^2,2));

% ------------------------------------------
% Calculate dual-control-volume cell areas
% ------------------------------------------
Area2 = zeros(n_verts,1);
x = zeros(4,2);
for iv=1:n_verts
    ic = v2c(iv,:);
    ic = ic(ic>0);
    Area2(iv) = sum(Area(ic)./c2nf(ic));
end

% ------------------------------------------
% Match up ghost nodes to corresponding interior node
% ------------------------------------------
n_bnd_verts = n_verts_total - n_verts;
n_bnd_edges = n_edges_total - n_edges;

% Ghost vertex to corresponding (mirror) interior vertex
gv2iv = zeros(n_bnd_verts,1);
% Ghost vertex to boundary condition of connected boundary edge
gv2bc = zeros(n_bnd_verts,1);
gv2be = zeros(n_bnd_verts,1);
for iv=1:n_bnd_verts
    % Find ghost cell (doesn't matter which) touching the ghost node
    icg = v2c(iv+n_verts,1);
    % Find the boundary edge of that ghost cell (index of bnde)
    bie = find(e2c(bnde,2)==icg);
    % Find connected vertex which lies on boundary
    vlist = v2v(iv+n_verts,:);
    ive = find(vlist<=n_verts);
    ive = ive(vlist(ive)>0);
    ive = v2v(iv+n_verts,ive);
    % Boundary condition of edge --> apply to ghost node
    gv2bc(iv) = bc(bie);
    gv2be(iv) = bnde(bie);
    
    % Now find the correct interior point
    % (look at cell across boundary edge and find non-boundary point to 
    %  left/right of the boundary point previously found)
    ie = bnde(bie);  % Global edge index
    ici = e2c(ie,1);
    k = find(c2v(ici,:)==ive);
    if isempty(k)
        error('sumthins wrong!');
    end
    
    if k==c2nf(ici)
        kp1 = 1;
    else
        kp1 = k+1;
    end
    if k==1
        km1 = c2nf(ici);
    else
        km1 = k-1;
    end
    
    if c2v(ici,kp1)==e2v(ie,1) || c2v(ici,kp1)==e2v(ie,2)
        % point is on bnde; wrong direction, grab other point
        gv2iv(iv) = c2v(ici,km1);
    else
        gv2iv(iv) = c2v(ici,kp1);
    end
end

ivi = [0,0];
ivg = [0,0];
for ie=1:n_bnd_edges
    ici = e2c(ie,1);
    icg = e2c(ie,2);
    iv1 = e2v(ie,1);
    iv2 = e2v(ie,2);
    
    ki=0;
    kg=0;
    for j=1:4
        iv0 = c2v(ici,j);
        if (iv0~=iv1) && (iv0~=iv2) % not on the boundary edge
            ki = ki + 1;
            ivi(ki) = iv0;
        end
        iv0 = c2v(icg,j);
        if (iv0~=iv1) && (iv0~=iv2) % not on the boundary edge
            kg = kg + 1;
            ivg(kg) = iv0;
        end
    end
end

save(meshFile,'xv','c2v','c2e','e2c','c2c','e2f','c2nf','e2v','c2n',...
    'unorm','dA','Area','XYC','XYE','dXE','bnde','inte','ie_sub_in','n_sub_in',...
    'ie_sub_out','n_sub_out','ie_sup_in','n_sup_in','ie_sup_out','n_sup_out',...
    'ie_slip_wall','n_slip_wall','ie_adiabatic_wall','n_adiabatic_wall',...
    'ie_isothermal_wall','n_isothermal_wall','n_cells','n_cells_total',...
    'n_edges','n_edges_total','bounds','bc','v2c','v2nc','v2ne','v2e','v2v','n_verts',...
    'n_bnd_edges','n_int_edges','max_nf','xmin','xmax','ymin','ymax',...
    'gv2iv','gv2bc','gv2be','dA2','Area2');

disp('  Mesh generation complete!');
disp(['  Mesh saved to ',meshFile,'.']);