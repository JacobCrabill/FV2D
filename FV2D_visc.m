%% AA 214C - Project 2 - Viscous Flat Plate
%  2nd-Order Finite Volume Scheme for the 2D Navier-Stokes Equations
%  Jacob Crabill - 5/27/14
clear all; close all; clc;
addpath('./data');
addpath('./mesh');

%% ------------------------------------------------------------------------
% Overall Simulation  Parameters
% -------------------------------------------------------------------------
% Plotting
plotMode = true;       % Plot continuously during run or not
plotFreq = 100;        % Frequecy at which to plot solution
resFreq = 100;         % Residual plotting frequency
plotVar = 6;           % Variable to plot (1-4: W, 5-8: phi, 9: M)
% Simulation Control
CFL = 2.4;             % CFL number
iterMax = 300000;       % Maximum iteration
limiter = 1;           % Slope limiter - Venkatakrishnan, Barth & Jespersen, or AJ
K = 1000;               % Venkatakrishnan limiter parameter (~.1-4? Lower = more limiting)
% Simulation Start & Input
restart = false;
meshFile    = 'ViscBump_Fine';  % File containing all geometry data to start from
restartFile = 'ViscBump_';
% Convergence
err_tol = 10^-5;       % Convergence tolerance (used on density)
resnorm = false;       % Use actual residual, or normalize to iteration 1?
% Data Output
outputFreq = 1000;    % Frequency of restart-file writing
dataFile = 'ViscBump_Fine';

%% ------------------------------------------------------------------------
% Load Restart File or Mesh File
% -------------------------------------------------------------------------
if restart
    load(restartFile);
else
    load(meshFile);
end

%% ------------------------------------------------------------------------
% Initial Conditions
% -------------------------------------------------------------------------
if ~restart
    T_inf = 300;                     % K         | Temperature
    p_inf = 5.8614*10^4;             % Pa        | Static Pressure
    theta = 0*pi/180;                % radians   | Flow angle
    M_inf = 0.908;                     % --        | Mach number
    gam = 1.4;                       % --        | Specific heat ratio
    R = 287.06;                      % J/kg-K    | Gas constant
    Cp = 1005;                       % J/kg-K    | Specific heat constant p
    Pr = .76;                        % --        | Prandtl number
end
rho_inf = p_inf/(R*T_inf);         % kg/m^3    | Density
a_inf = sqrt(gam*R*T_inf);         % m/s       | Speed of sound
n_inf = [cos(theta);sin(theta)];   % --        | Flow direction
u_inf = M_inf*a_inf*n_inf(1);      % m/s       | x-direction velocity
v_inf = M_inf*a_inf*n_inf(2);      % m/s       | y-direction velocity
E_inf = p_inf/(gam-1) + rho_inf*(M_inf*a_inf)^2/2; % J/m^3 = kg/m-s^2 | Energy
Tt_inf = T_inf*(1+(gam-1)/2*M_inf^2);  % K   | Total temperature
pt_inf = p_inf*(1+(gam-1)/2*M_inf^2)^(gam/(gam-1));  % Pa  | Total pressure
mu_ref = 1.716E-5; % Reference Viscosity
T_s = 110.56;      % Sutherland's Law Reference Temperature
T_ref = 273.11;    % Reference Temperature

%% ------------------------------------------------------------------------
% Solution / Flux Storage Initialization
% -------------------------------------------------------------------------
% Solution Variables
if ~restart
    W = zeros(n_cells_total,4); % Conserved variables
end
phi = zeros(n_cells_total,4);   % Primitive variables
WL = zeros(4,1);                % Temporary left-state for boundary edge
WR = zeros(4,1);                % Temporary right-state for boundary edge
WeL = zeros(n_edges,4);         % Extrapolated left state at all edges
WeR = zeros(n_edges,4);         % Extrapolated right state at all edges
WLi = zeros(n_int_edges,4);     % Left state for interface flux calculation
WRi = zeros(n_int_edges,4);     % Right state for interface flux calculation
WE = zeros(n_cells_total,max_nf,4);   % Conserved variables at faces of each cell
phiL = zeros(n_edges,4);        % Left state at edges (primitive variables)
phiR = zeros(n_edges,4);        % Right state at edges (primitive variables)
phiE = zeros(n_cells_total,max_nf,4); % Primitives at each face of each cell
phiN = zeros(n_cells_total,max_nf,4); % Primitives at each cell's neighbors

% Gradients
dWdX = zeros(4,2,n_cells);      % n_fields x n_dims x n_cells
dphidX = zeros(4,2,n_cells_total);    % n_fields x n_dims x n_cells
dphi = zeros(4,max_nf,n_cells);
% Limiter Variables
PHI = zeros(4,max_nf,n_cells_total);
PHIj = zeros(4,n_cells_total);
PHImin = zeros(4,n_cells_total);

% --- Viscous Terms ---
% Solution Variables
tau = zeros(n_cells_total,6);   % viscous stress terms (plus temp.)
tauL = zeros(n_edges,6);        % Left state at edges (primitive variables)
tauR = zeros(n_edges,6);        % Right state at edges (primitive variables)
tauE = zeros(n_cells_total,max_nf,6); % Primitives at each face of each cell
tauN = zeros(n_cells_total,max_nf,6); % Primitives at each cell's neighbors
% Gradients
dtaudX = zeros(6,2,n_cells);    % n_fields x n_dims x n_cells
% Limiter Variables
TAU = zeros(6,max_nf,n_cells_total);
TAUj = zeros(6,n_cells_total);
TAUmin = zeros(6,n_cells_total);

% Other
fnorm = zeros(2,1);
RES = zeros(n_cells,4);
Fn = zeros(n_edges,4);  % F (dot) n
if ~restart
    res = ones(iterMax,4);      % Solution residual time history
else
    res = [res;ones(iterMax-iter,4)];
end

% RK time-stepping variables
W0 = W;                    % Solution at step 0    
WN = zeros(n_cells,4);     % Partially-assembled W^(n+1)
a_rk = [1/6,1/3,1/3,1/6];  % RK4 "a" constants
b_rk = [0, 1/2, 1/2, 1];   % RK4 "b" constants

% ------------------------------------------------------------------------
% Initialize Solution
% ------------------------------------------------------------------------
if ~restart
    W(:,1) = rho_inf;
    W(:,2) = rho_inf*u_inf;
    W(:,3) = rho_inf*v_inf;
    W(:,4) = E_inf;
end

% ------------------------------------------------------------------------
% Store Least-Squares Gradient-Estimate Matrices
% ------------------------------------------------------------------------
% Setup QR factorization of gradient estimate equations
max_nc0 = 15;
DX = zeros(max_nc0,2,n_cells);
Q1 = zeros(max_nc0,max_nc0,n_cells); % actually Q' for efficiency
R1_qr = zeros(max_nc0,2,n_cells);
c2ac1 = zeros(n_cells,max_nc0);
c2nac = zeros(n_cells,1);

max_nc = 0;
for ic=1:n_cells
    ic2 = unique(v2c(c2v(ic,:),:));
    ic2 = ic2(ic2>0);
    ic2 = ic2(ic2~=ic);
    nc = length(ic2);
    max_nc = max(max_nc,nc);
    dX = zeros(nc,2);
    for j=1:nc
        dX(j,:) = XYC(ic2(j),:) - XYC(ic,:);
    end
    [q,r] = qr(dX);
    Q1(1:nc,1:nc,ic) = q';
    R1_qr(1:nc,:,ic) = r;
    c2ac1(ic,1:nc) = ic2;
    c2nac(ic) = nc;
end

Q = zeros(max_nc,max_nc,n_cells); % actually Q' for efficiency
R_qr = zeros(max_nc,2,n_cells);
dXinv = zeros(2,max_nc,n_cells);
c2ac = zeros(n_cells,max_nc);
if max_nc < max_nc0
    for ic=1:n_cells
        Q(:,:,ic) = Q1(1:max_nc,1:max_nc,ic);
        R_qr(:,:,ic) = R1_qr(1:max_nc,:,ic);
        dXinv(:,:,ic) = R_qr(:,:,ic)\Q(:,:,ic);
        c2ac(ic,:) = c2ac1(ic,1:max_nc);
    end
end

% ------------------------------------------------------------------------
% Calculate Misc. Constants
% ------------------------------------------------------------------------
astar2 = 2*gam*R/(gam+1)*Tt_inf;
tanTheta = tan(theta);
tanTheta2 = tanTheta^2;
dpdu_1 = (gam-1)/(gam+1)*(tanTheta2+1)/astar2;
dpdu_2 = -2/(gam+1)*(1+tanTheta2)/astar2*gam*pt_inf;
delta_x = sqrt(4/sqrt(3)*sum(Area(1:n_cells))/n_cells);
eps2 = (K*delta_x).^3;
n_verts_total = length(xv);
n_edges_total = length(e2v);
n_cells_total = length(c2v);

% ------------------------------------------------------------------------
% Plotting & Start Time
% ------------------------------------------------------------------------
if plotMode
    h1 = figure(1);
    n_cells_plot = n_cells;
    ic = (1:n_cells_plot)';
    x = reshape(xv(c2v(ic,:),1),n_cells_plot,max_nf)';
    y = reshape(xv(c2v(ic,:),2),n_cells_plot,max_nf)';
    Vars = cell(9);
    Vars{1} = 'Density';
    Vars{2} = 'X-Momentum';
    Vars{3} = 'Y-Momentum';
    Vars{4} = 'Energy';
    Vars{5} = 'Density';
    Vars{6} = 'X-Velocity';
    Vars{7} = 'Y-Velocity';
    Vars{8} = 'Pressure';
    Vars{9} = 'Mach Number';
end
if resFreq >0; h2 = figure(2); end

%% ------------------------------------------------------------------------
% Solve
% -------------------------------------------------------------------------
%tic;
if ~restart; iter = 0; end
while iter < iterMax
    % RK4 Time-Stepping
    W0 = W;
    WN = W0(1:n_cells,:);
    Fn = 0*Fn;
    for rk_i = 1:4
        % Calculate primitive variables
        phi(:,1) = W(:,1);                                            % rho
        phi(:,2) = W(:,2)./W(:,1);                                      % u
        phi(:,3) = W(:,3)./W(:,1);                                      % v
        phi(:,4) = (gam-1)*(W(:,4)-0.5*(W(:,2).^2+W(:,3).^2)./W(:,1));  % p
                
        % ================ INVISCID BOUNDARY CONDITIONS ================= %
        % --------------- SUBSONIC INLET --------------- %
        if n_sub_in > 0
            ie = ie_sub_in;
            ici = e2c(ie,1);
            icg = e2c(ie,2);
            fnorm = unorm(ie,:);
            da = dA(ie);
            rho = phi(ici,1);
            u = phi(ici,2);
            v = phi(ici,3);
            p = phi(ici,4);
            E = p/(gam-1) + 0.5*rho.*(u.^2+v.^2);
            a = sqrt(gam*p./rho);

            % ---------------------------------------------
            % Get cell "(2,j)" from cell "(1,j)"
            % NOTE: This will only work for quads!
            wallFace = e2f(ie,1); % local face ID of interior cell
            wallFace = wallFace - 2; % face opposite wall
            wallFace(wallFace<0) = wallFace(wallFace<0) + 4;
            ic2 = zeros(size(ici));
            for i=1:n_sub_in
                ic2(i) = c2c(ici(i),wallFace(i));
            end
            dp = phi(ic2,4) - phi(ici,4);
            du = phi(ic2,2) - phi(ici,2);
            % ---------------------------------------------
            
            dt = CFL*dA(ie)./(sqrt(u.^2+v.^2)+a);
            dX = XYC(ici,:)-XYC(icg,:); % 'right-to-left' is correct: see MacCormack
            dx = sqrt(sum(dX.^2,2));
            
            lambda = (u-a).*dt./dx;
            A = (dp-rho.*a.*du).*lambda./(lambda-1);
            dpdu = dpdu_2*u.*(1-dpdu_1*u.^2).^(1/(gam-1));

            deltau = A./(dpdu - rho.*a);
            un = u+deltau;
            vn = tanTheta*un;
            stag_ratio = (1-dpdu_1*u.^2);
            Tn = Tt_inf*stag_ratio;
            pn = pt_inf*stag_ratio.^(gam/(gam-1));            
            %rhon = pn/(R*Tt_inf);  % ??? was this correct? doesn't seem so
            rhon = pn./(R*Tn);
            
            En = pn/(gam-1) + 0.5*rhon.*(un.^2+vn.^2);
            
            WL = W(ici,:);
            WR = zeros(length(ie),4);
            WR(:,1) = rhon;
            WR(:,2) = rhon.*un;
            WR(:,3) = rhon.*vn;
            WR(:,4) = En;
            
            % Right state goes to ghost cell
            W(icg,:) = WR;
            phi(icg,:) = [rhon,un,vn,pn];

            clear FnIE;
            FnIE = convective_flux(WL,WR,fnorm,da);
            Fn(ie,:) = FnIE;
        end
        % --------------- SUBSONIC OUTLET --------------- %
        if n_sub_out > 0
            ie = ie_sub_out;
            ici = e2c(ie,1);
            icg = e2c(ie,2);
            fnorm = unorm(ie,:);
            da = dA(ie);
            rho = phi(ici,1);
            u = phi(ici,2);
            v = phi(ici,3);
            p = phi(ici,4);     
            E = p/(gam-1) + 0.5*rho.*(u.^2+v.^2);   
            a = sqrt(gam*p./rho);

            % ---------------------------------------------
            % Get cell "(2,j)" from cell "(1,j)"
            % NOTE: This will only work for quads!
            wallFace = e2f(ie,1); % local face ID of interior cell
            wallFace = wallFace - 2; % face opposite wall
            wallFace(wallFace<=0) = wallFace(wallFace<=0) + 4;
            ic2 = zeros(size(ici));
            for i=1:n_sub_out
                ic2(i) = c2c(ici(i),wallFace(i));
            end
            drho = phi(ic2,1) - phi(ici,1);
            dp = phi(ic2,4) - phi(ici,4);
            du = phi(ic2,2) - phi(ici,2);
            dv = phi(ic2,3) - phi(ici,3);
            % ---------------------------------------------
           
            dt = CFL*dA(ie)./(sqrt(u.^2+v.^2)+a);
            dX = XYC(ici,:)-XYC(icg,:); % 'right-to-left' is correct: see MacCormack
            dx = sqrt(sum(dX.^2,2));
            
            lambda1 = u.*dt./dx;
            lambda2 = (u+a).*dt./dx;
            lambda4 = (u-a).*dt./dx;

            R1 = -lambda1./(lambda1+1).*(drho-dp./a.^2);
            R2 = -lambda2./(lambda2+1).*(dp+rho.*a.*du);
            R3 = -lambda1./(lambda1+1).*dv;
            R4 = -lambda4./(lambda4+1).*(dp-rho.*a.*du);

            deltaP = zeros(length(ie),1);
            ip = find(abs(u./a)>1);  % M > 1; must specify all quantities
            deltaP(ip) = (R2(ip)+R4(ip))/2;

            deltarho = R1 + deltaP./a.^2;
            deltau = (R2-deltaP)./(rho.*a);

            pn = p + deltaP;
            rhon = rho + deltarho;
            un = u + deltau;
            vn = v + R3;
            En = pn/(gam-1) + 0.5*rhon.*(un.^2+vn.^2);

            WL = W(ici,:);
            WR = zeros(length(ie),4);
            WR(:,1) = rhon;
            WR(:,2) = rhon.*un;
            WR(:,3) = rhon.*vn;
            WR(:,4) = En;

            % Right state goes to ghost cell
            W(icg,:) = WR;
            phi(icg,:) = [rhon,un,vn,pn];

            clear FnIE;
            FnIE = convective_flux(WL,WR,fnorm,da);
            Fn(ie,:) = FnIE;
        end
        % --------------- SUPERSONIC INLET --------------- %
        if n_sup_in > 0
            ie = ie_sup_in;
            fnorm = unorm(ie,:);
            
            % Euler flux
            rhou_n = rho_inf*[u_inf,v_inf]*fnorm';
            Fn(ie,1) = rhou_n;
            Fn(ie,2:3) = rhou_n'*[u_inf,v_inf] + p_inf*fnorm;
            Fn(ie,4) = rhou_n*(E_inf + p_inf)/rho_inf;
            Fn(ie,:) = Fn(ie,:).*repmat(dA(ie),[1,4]);
        end
        % --------------- SUPERSONIC OUTLET --------------- %
        if n_sup_out > 0
            ie = ie_sup_out;
            ici = e2c(ie,1);
            icg = e2c(ie,2);
            fnorm = unorm(ie,:);

            WL = W(ici,:);
            rhou_n = WL(:,2).*fnorm(:,1) + WL(:,3).*fnorm(:,2);

            % Update ghost cells
            W(icg,:) = WL;
            phi(icg,:) = phi(ici,:);
            tau(icg,:) = tau(ici,:);
            
            % Euler flux
            Fn(ie,1) = rhou_n;
            Fn(ie,2) = rhou_n.*WL(:,2)./WL(:,1) + phi(ici,4).*fnorm(:,1);
            Fn(ie,3) = rhou_n.*WL(:,3)./WL(:,1) + phi(ici,4).*fnorm(:,2);
            Fn(ie,4) = rhou_n.*(WL(:,4) + phi(ici,4))./WL(:,1);
            Fn(ie,:) = Fn(ie,:).*repmat(dA(ie),[1,4]);
        end
        % --------------- SLIP WALL --------------- %
        if n_slip_wall > 0
            ie = ie_slip_wall;
            ici = e2c(ie,1);
            icg = e2c(ie,2);

            pwall = phi(ici,4);
            
            if find(pwall<0,1)
                disp('WARNING: Negative pressure encountered.');
                p_tol = min(pwall(pwall>0));
                pwall(pwall<0) = p_tol; %1E-6;
            end

            % Reflect velocity about edge (U_r = U - 2*(U dot nhat))
            ug = phi(ici,2).*(1-2*abs(unorm(ie,1)));
            vg = phi(ici,3).*(1-2*abs(unorm(ie,2)));

            % Assign quantities to ghost cell
            phi(icg,1) = phi(ici,1);
            phi(icg,2) = ug;
            phi(icg,3) = vg;
            phi(icg,4) = pwall;
            W(icg,1) = W(ici,1);
            W(icg,2) = W(ici,1).*ug;
            W(icg,3) = W(ici,1).*vg;
            W(icg,4) = pwall/(gam-1) + 0.5*phi(ici,1).*(ug.^2+vg.^2);

            % Extrapolate pressure to wall [aka use reconstructed left state]        
            Fn(ie,1) = 0;
            Fn(ie,2) = pwall.*unorm(ie,1);
            Fn(ie,3) = pwall.*unorm(ie,2);
            Fn(ie,4) = 0;
            Fn(ie,:) = Fn(ie,:).*repmat(dA(ie),[1,4]);
        end
        % --------------- ADIABATIC WALL --------------- %
        if n_adiabatic_wall > 0
            ie = ie_adiabatic_wall;
            ici = e2c(ie,1);
            icg = e2c(ie,2);
            
            pwall = phi(ici,4);
            rho = phi(ici,1);           

            % Reflect velocity about edge (U_r = U - 2*(U dot nhat))
            ui =  phi(ici,2);
            vi =  phi(ici,3);
            ug = -phi(ici,2);
            vg = -phi(ici,3);
                       
            uwall = zeros(size(ui));
            phi_l = [rho,uwall,uwall,pwall];
            phi_r = [rho,uwall,uwall,pwall];

            % Assign quantities to ghost cell
            phi(icg,1) = rho;
            phi(icg,2) = ug;
            phi(icg,3) = vg;
            phi(icg,4) = pwall;
            W(icg,1) = rho;
            W(icg,2) = rho.*ug;
            W(icg,3) = rho.*vg;
            W(icg,4) = pwall/(gam-1) + 0.5*rho.*(ug.^2+vg.^2);

            % Extrapolate pressure to wall [aka use reconstructed left state]        
            Fn(ie,1) = 0;
            Fn(ie,2) = pwall.*unorm(ie,1);
            Fn(ie,3) = pwall.*unorm(ie,2);
            Fn(ie,4) = 0;
            Fn(ie,:) = Fn(ie,:).*repmat(dA(ie),[1,4]);
        end
        
        % ==================== GRADIENT CALCULATION ===================== %
        % Use a least-squares fit through surrounding cells
        dphidX = calc_grad_dxinv(phi,c2ac,c2nac,dXinv);
        
        % ------------------- Apply B.C.'s to gradients -------------------
%         icg = ((n_cells+1):n_cells_total)';
%         ici = c2c(icg,1);
%         bcg = bc(c2e(icg,1));
%         
%         % Inlet/Outlet - Copy
%         icg1 = icg(bcg<5);
%         ici1 = ici(bcg<5);
%         dphidX(:,:,icg1) = dphidX(:,:,ici1);
%         
%         % NOTE: ghost cells ALWAYS on 'right' of all boundary edges
%         % Slip Wall - Mirror velocities about wall
%         if n_slip_wall > 0
%             icg1 = e2c(ie_slip_wall,2);
%             ici1 = e2c(ie_slip_wall,1);
%             % du/dx and dv/dy are the same
%             dphidX(2,1,icg1) = dphidX(2,1,ici1);
%             dphidX(3,2,icg1) = dphidX(3,2,ici1);
%             % du/dy and dv/dx are mirrored
%             dphidX(2,2,icg1) = -dphidX(2,2,ici1);
%             dphidX(3,1,icg1) = -dphidX(3,1,ici1);
%             % drho/dx and dp/dx are same
%             dphidX([1,4],1,icg1) = dphidX([1,4],1,ici1);
%             % drho/dy and dp/dy are mirrored
%             dphidX([1,4],2,icg1) = -dphidX([1,4],2,ici1);
%         end
%         % No-Slip Wall - Mirror and flip gradients about wall
%         % Do coordinate transformation at edges from x,y to  normal, 
%         % tangential coordinates
%         if n_adiabatic_wall > 0 || n_isothermal_wall > 0
%             icg1 = [e2c(ie_adiabatic_wall,2);e2c(ie_isothermal_wall,2)];
%             ici1 = [e2c(ie_adiabatic_wall,1);e2c(ie_isothermal_wall,1)];
%             
%             ie = [ie_adiabatic_wall,ie_isothermal_wall];
%             nhat = unorm(ie,:);
%             blah = ones(length(nhat),2)/sqrt(2);
%             that = blah - repmat(sum(nhat.*blah,2),[1,2]).*nhat;
%             that = that./repmat(sqrt(sum(that.^2,2)),[1,2]);
%             dundu = nhat(:,1); dundv = nhat(:,2);
%             dutdu = that(:,1); dutdv = that(:,2);
%             dxdn = 1./nhat(:,1); dydn = 1./nhat(:,2);
%             dxdt = 1./that(:,1); dydt = 1./that(:,2);
%             % dudn = dudx*dxdn + dudy*dydn, etc.
%             dudn = squeeze(dphidX(2,1,ici1)).*dxdn + squeeze(dphidX(2,2,ici1)).*dydn;
%             dvdn = squeeze(dphidX(3,1,ici1)).*dxdn + squeeze(dphidX(3,2,ici1)).*dydn;
%             dudt = squeeze(dphidX(2,1,ici1)).*dxdt + squeeze(dphidX(2,2,ici1)).*dydt;
%             dvdt = squeeze(dphidX(3,1,ici1)).*dxdt + squeeze(dphidX(3,2,ici1)).*dydt;
%             drhodn = squeeze(dphidX(1,1,ici1)).*dxdn + squeeze(dphidX(1,2,ici1)).*dydn;
%             dpdn   = squeeze(dphidX(4,1,ici1)).*dxdn + squeeze(dphidX(4,2,ici1)).*dydn;
%             drhodt = squeeze(dphidX(1,1,ici1)).*dxdt + squeeze(dphidX(1,2,ici1)).*dydt;
%             dpdt   = squeeze(dphidX(4,1,ici1)).*dxdt + squeeze(dphidX(4,2,ici1)).*dydt;
%             % dun/dn = dun/du * du/dn + dun/dv * dv/dn, etc.
%             % dut/dn, dun/dn are same
%             % dut/dt, dun/dt are mirrored
%             dundn = dundu.*dudn + dundv.*dvdn;
%             dundt =-dundu.*dudt - dundv.*dvdt;
%             dutdn = dutdu.*dudn + dutdv.*dvdn;
%             dutdt =-dutdu.*dudt - dutdv.*dvdt;
%             
%             % Transform un, ut to x,y
%             dundx = dundn./dxdn + dundt./dxdt;
%             dundy = dundn./dydn + dundt./dydt;
%             dutdx = dutdn./dxdn + dutdt./dxdt;
%             dutdy = dutdn./dydn + dutdt./dydt;
%             
%             % Transform all derivatives back to x,y
%             dphidX(2,1,icg1) = (1./dundu).*dundx + (1./dutdu).*dutdx; % dudx
%             dphidX(2,2,icg1) = (1./dundu).*dundy + (1./dutdu).*dutdy; % dudy
%             dphidX(3,1,icg1) = (1./dundv).*dundx + (1./dutdv).*dutdx; % dvdx
%             dphidX(3,2,icg1) = (1./dundv).*dundy + (1./dutdv).*dutdy; % dvdy
%             % drho/dt, dp/dt are same
%             % drho/dn, dp/dn are mirrored
%             dphidX(1,1,icg1) = -drhodn./dxdn + drhodt./dxdt; % drhodx
%             dphidX(1,2,icg1) = -drhodn./dydn + drhodt./dydt; % drhody
%             dphidX(4,1,icg1) = -dpdn./dxdn + dpdt./dxdt;     % dpdx
%             dphidX(4,2,icg1) = -dpdn./dydn + dpdt./dydt;     % dpdy
%         end
        size(dphidX)
        temp = apply_bcs_grad_phi(dphidX,c2c,c2e,bc,unorm,n_cells);
        dphidX = temp;
        
        % ====================== GRADIENT LIMITING ====================== %
        % Reconstruct solution: extrapolate from cell centers to edge centers
        dphi = calc_dU(dphidX,dXE,n_cells_total);
        
        ic = (1:n_cells_total)';
        for j = 1:max_nf
            phiE(ic,j,:) = phi(ic,:) + squeeze(dphi(:,j,ic))';
        end

        % Get min & max of neighboring cells
        for j = 1:max_nf
            ic2 = c2c(ic,j);
            ic2 = ic2(ic2>0);
            phiN(ic(ic2>0),j,:) = phi(ic2,:);
        end
        phiMax = squeeze(max(phiN,[],2));
        phiMin = squeeze(min(phiN,[],2));

        % Apply limiting to the primitive variables
        % Venkatakrishnan's Limiter
        ic = (1:n_cells_total)';
        D1_max = phiMax - phi(ic,:);
        D1_min = phiMin - phi(ic,:);

        for j = 1:max_nf
            for k = 1:4
                dphijk = squeeze(dphi(k,j,:));

                ic1 = find(dphijk>0);
                D1maxk = D1_max(ic1,k);
                D2 = 0.5*dphijk(ic1);
                PHI(k,j,ic1) = ((D1maxk.^2+eps2).*D2 + 2*D2.^2.*D1maxk)./(D2.*(D1maxk.^2 + 2*D2.^2 + D1maxk.*D2 + eps2));

                ic2 = find(dphijk<0);
                D1mink = D1_min(ic2,k);
                D2 = 0.5*dphijk(ic2);
                PHI(k,j,ic2) = ((D1mink.^2+eps2).*D2 + 2*D2.^2.*D1mink)./(D2.*(D1mink.^2 + 2*D2.^2 + D1mink.*D2 + eps2));

                ic3 = find(dphijk==0);
                PHI(k,j,ic3) = 1;
            end
        end
        
        % Get the final limiting function & apply to del(phi)
        ic = (1:n_cells_total)';
        PHImin = min(PHI,[],2);
%         ic1 = e2c(ie_adiabatic_wall,:);
%         PHImin(2:3,1,ic1) = 1;
        for j = 1:2
            dphidX(:,j,ic) = dphidX(:,j,ic) .* PHImin;
        end
        
%         ic1 = e2c(ie_adiabatic_wall,1);
%         ic2 = e2c(ie_adiabatic_wall,2);
%         dphidX(2,2,ic1) = phi(ic1,2)./(XYC(ic1,2)-XYE(ie_adiabatic_wall,2));
%         dphidX(3,2,ic1) = phi(ic1,3)./(XYC(ic1,2)-XYE(ie_adiabatic_wall,2));
%         dphidX(2,2,ic2) = phi(ic2,3)./(XYC(ic2,2)-XYE(ie_adiabatic_wall,2));
%         dphidX(3,2,ic2) = phi(ic2,4)./(XYC(ic2,2)-XYE(ie_adiabatic_wall,2));
        
        % ======================== VISCOUS TERMS ======================== %
        % Calculate gradients of viscous terms for reconstruction
        ic = (1:n_cells_total)';
        tau(ic,1:2) = squeeze(dphidX(2,1:2,ic))';  % dudx, dudy
        tau(ic,3:4) = squeeze(dphidX(3,1:2,ic))';  % dvdx, dvdy
        % dT/dx = R/p drho/dx - r*rho/p^2 dp/dx    [Temperature]
        tau(ic,5) = 1./(R*phi(ic,1)) .* squeeze(dphidX(4,1,ic)) -...
                   phi(ic,4)./(R*phi(ic,1).^2) .* squeeze(dphidX(1,1,ic));
        % dT/dy = R/p drho/dy - r*rho/p^2 dp/dy    [Temperature]
        tau(ic,6) = 1./(R*phi(ic,1)) .* squeeze(dphidX(4,2,ic)) -...
                   phi(ic,4)./(R*phi(ic,1).^2) .* squeeze(dphidX(1,2,ic));
               
        % ===================== SECOND DERIVATIVES ====================== %
        dtaudX = calc_grad_tau(tau,v2c,v2nc,xv,e2v,v2ne,c2v,c2nf,gv2iv,gv2bc);
        % Use a least-squares fit through surrounding cells             
%         dtaudX = zeros(6,2,n_cells_total);
%         for ic=1:n_cells
%             ic2 = unique(v2c(c2v(ic,:),:));
%             ic2 = ic2(ic2>0);
%             nc = length(ic2);
%             dtau = zeros(nc,6);
%             dX = zeros(nc,2);
%             for j=1:nc
%                 dtau(j,:) = tau(ic2(j),:) - tau(ic,:);
%                 dX(j,:) = XYC(ic2(j),:) - XYC(ic,:);
%             end
%             dtaudX(:,:,ic) = (dX\dtau)';
%             dtaudX(:,:,ic) = (R_qr(1:nc,:,ic)\Q(1:nc,1:nc,ic)*dtau)';
%             dtaudX(:,:,ic) = (dXinv(:,1:nc,ic)*dtau)';
%         end

        % Mirror second derivatives to ghost cells
        icg = (n_cells+1):n_cells_total;
        ici = c2c(icg,1);
        dtaudX(:,:,icg) = dtaudX(:,:,ici);
        
        % ================== VISCOUS GRADIENT LIMITING ================== %
        % Reconstruct solution: extrapolate from cell centers to edge centers
        dtau = calc_dU(dtaudX,dXE,n_cells_total);
                        
        % =================== SOLUTION RECONSTRUCTION =================== %
        % Reconstruct solution: extrapolate from cell centers to edge centers
        % Apply the limiter to the previously calculated gradients
        % get the new reconstructed phi and tau vectors
        ic = (1:n_cells_total)';
        for j = 1:max_nf
            phiE(ic,j,:) = phi(ic,:) + squeeze(PHImin.*dphi(:,j,ic))';
%             phiE(ic,j,:) = phi(ic,:) + squeeze(dphi(:,j,ic))';
            tauE(ic,j,:) = tau(ic,:) + squeeze(dtau(:,j,ic))';
        end

        % Calculate new conservative variables at edge of each element
        WE(:,:,1) = phiE(:,:,1);
        WE(:,:,2) = phiE(:,:,1).*phiE(:,:,2);
        WE(:,:,3) = phiE(:,:,1).*phiE(:,:,3);
        WE(:,:,4) = phiE(:,:,4)/(gam-1) + 0.5*phiE(:,:,1).*(phiE(:,:,2).^2+phiE(:,:,3).^2);

        % =================== LEFT & RIGHT EDGE STATES ================== %
        % Reorganize reconstructed solution into L & R states for each edge
        % Need to match local face of each cell to global edge
        % WL(ie,:) = WE(icL,ifL,:);
        % WR(ie,:) = WE(icR,ifR,:);
        ie = (1:n_edges)';
        icL = e2c(ie,1);
        icR = e2c(ie,2);

        ifL = e2f(ie,1);
        ifR = e2f(ie,2);

        for j = 1:max_nf
            IEL = find(ifL==j);
            ICL = icL(IEL);
            WeL(IEL,:) = WE(ICL,j,:);
            phiL(ie(IEL),:) = phiE(ICL,j,:);
            tauL(ie(IEL),:) = tauE(ICL,j,:);

            IER = find(ifR==j);
            ICR = icR(IER);
            WeR(IER,:) = WE(ICR,j,:);
            phiR(ie(IER),:) = phiE(ICR,j,:);
            tauR(ie(IER),:) = tauE(ICR,j,:);
        end        
        WLi = WeL(inte,:);
        WRi = WeR(inte,:);
        
        % ======================= CONVECTIVE FLUX ======================= %
        % Convective portion of flux already specified at the boundaries,
        % so only needed at interior interfaces
        fnorm = unorm(inte,:);
        da = dA(inte);
        
        phiLi = phiL(inte,:);
        phiRi = phiR(inte,:);
                
        clear FnIE;
        FnIE = Roe_Flux_2(WLi,WRi,fnorm,da);
        Fn(inte,:) = FnIE;
        
        % ======================== VISCOUS FLUX ========================= %
        % Viscous flux not yet calculated anywhere, so needed at all
        % interfaces
        % err, try to enforce a no-slip wall:
        %phiL(ie_adiabatic_wall,2) = 0;  phiR(ie_adiabatic_wall,2) = 0;
        %phiL(ie_adiabatic_wall,3) = 0;  phiR(ie_adiabatic_wall,3) = 0;
        
        da = dA(1:n_edges);
        clear FnVis;
        FnVis = viscous_flux(tauL,tauR,phiL,phiR,unorm,da,Cp,Pr,R);
        Fn = Fn + FnVis;
        
        % ====================== TIME ADVANCEMENT ======================= %
        % Estimate convective time step limit
        ic =(1:n_cells)';
        a = sqrt(gam*phi(ic,4)./phi(ic,1));
        magU = sqrt(phi(ic,2).^2+phi(ic,3).^2);
        DTc = CFL*min(dA(c2e(ic,:)),[],2)./(magU+a);
        
        % Estimate viscous time step limit
        T = phi(ic,4)./(phi(ic,1)*R);
        nu = mu_ref*(T/T_ref).^(3/2)*(T_ref+T_s)./(T+T_s).*phi(ic,1);
        %DTv = CFL*min(dA(c2e(ic,:)).^2,[],2)./nu;
        %DTv = CFL*min(dA(c2e(ic,:))/2,[],2)./nu;
        
        % Get minimum time step
        %DT = min(DTc,DTv);
        DT = repmat(DTc,[1,max_nf,4]);
        
        % Gather flux contributions for each element
        FN = reshape(Fn(c2e(ic,:),:),n_cells,max_nf,4); % [cell,face,field]
        C2N = repmat(c2n(ic,:),[1,1,4]);
        j = 1:c2nf(ic);
        
        % Calculate the residual
        RES = reshape(sum(DT.*C2N.*FN,2),n_cells,4)./repmat(Area(ic),[1,4]);
        
        % RK time-stepping - W is technically "K_(i+1)"; WN adds up to W^(n+1)
        if rk_i<4
            W(ic,:) = W0(ic,:) - b_rk(rk_i+1)*RES;  % K_{i+1}
        end
        WN = WN - a_rk(rk_i)*RES;  % Running total of W^{n+1}
    end
    W(ic,:) = WN;
    
    iter = iter + 1;
    
    % ======================== SOLUTION RESIDUAL ======================== %
    if (resnorm)
        res(iter,:) = sum(abs(RES),1)/n_cells./res(1,:);
    else
        res(iter,:) = sum(abs(RES),1)/n_cells;
    end
    if mod(iter,resFreq)==0 && resFreq > 0
        set(0,'CurrentFigure',h2);
        semilogy(res,'linewidth',2); xlabel('Iteration'); ylabel('residual');
        grid on; xlim([0 iter]); ylim([min(min(res(2:iterMax,:)))/10,1]);
        legend('\rho','\rhou','\rhov','E');
        title(['Residual, iteration ',num2str(iter)]); drawnow;
    end
        
    % ======================== SOLUTION PLOTTING ======================== %
    if mod(iter,plotFreq)==0
        %toc
        if plotMode
            ic = (1:n_cells_plot)';
            set(0,'CurrentFigure',h1);
            if plotVar <= 4       % conserved variables
                c = W(ic,plotVar)'; 
            elseif plotVar <= 8   % primitive variables
                c = phi(ic,plotVar-4)'; 
            elseif plotVar == 9   % Mach number
                asq = gam*phi(ic,4)./W(ic,1);
                M = sqrt(sum(phi(ic,2:3).^2,2)./asq);
                c = M(ic)'; 
            end
            patch(x,y,c,'EdgeColor','none');
            grid on,xlabel('x'),ylabel('y'),colorbar;
            %xlim([xmin,xmax]),ylim([ymin,ymax]);
            xlim([-.1,.2]),ylim([0,.06]);
            var = Vars{plotVar};
            title([var,', Iteration ',num2str(iter)]), drawnow;
        end
        disp(['Iteration: ',num2str(iter)]);
    end
      
    % ========================== RESTART FILE =========================== %
    if mod(iter,outputFreq)==0
        fileName = ['./data/',dataFile,'_',num2str(iter)];
        disp(['Saving data to ',fileName]);
        save(fileName,'c2e','e2c','e2v','c2v','c2nf','c2n','c2c','xv',...
          'e2f','Area','dA','XYC','XYE','max_nf','bnde','inte','bounds',...
          'bc','unorm','ie_sub_in','ie_sub_out','ie_sup_in','ie_sup_out',...
          'ie_slip_wall','ie_isothermal_wall','n_sub_in','n_sub_out',...
          'ie_adiabatic_wall','n_adiabatic_wall','n_sup_in','n_sup_out',...
          'n_slip_wall','n_isothermal_wall','v2c','v2nc','v2ne',...
          'n_cells','n_cells_total','n_edges','n_verts','xmax','xmin',...
          'ymax','ymin','M_inf','p_inf','theta','T_inf','gam','R','Cp',...
          'Pr','W','phi','phiL','phiR','phiE','PHImin','PHI','dXE',...
          'dphidX','tau','dtaudX','tauL','tauR','tauE','TAUmin','TAU',...
          'Fn','dWdX','res','iter','n_int_edges','n_bnd_edges','gv2iv',...
          'gv2bc','gv2be','dA2','Area2')
    end
    
    if find(isnan(W(:,1)))
        error(['NaN! At iteration ',num2str(iter)]);
    end
end
