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
CFL = .5;             % CFL number
iterMax = 300000;       % Maximum iteration
limiter = 1;           % Slope limiter - Venkatakrishnan, Barth & Jespersen, or AJ
K = 1000;               % Venkatakrishnan limiter parameter (~.1-4? Lower = more limiting)
% Simulation Start & Input
restart = false;
meshFile    = 'testBox';  % File containing all geometry data to start from
restartFile = 'ViscBump_MC_K1000_12000';
% Convergence
err_tol = 10^-5;       % Convergence tolerance (used on density)
resnorm = false;       % Use actual residual, or normalize to iteration 1?
% Data Output
outputFreq = 1000;    % Frequency of restart-file writing
dataFile = 'testBox';

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

% ==================== GRADIENT CALCULATION ===================== %
% Use a least-squares fit through surrounding cells
phi(:,1) = XYC(:,1);
phi(:,2) = sum(XYC,2);
phi(:,3) = 2*XYC(:,1)+XYC(:,2);

tic;
% The mex'd way
dphidX1 = calc_grad_dxinv(phi,c2ac,c2nac,dXinv);
toc
tic;
% The slow way
dphidX = zeros(4,2,n_cells_total);        
for ic=1:n_cells
    ic2 = unique(v2c(c2v(ic,:),:));
    ic2 = ic2(ic2>0);
    ic2 = ic2(ic2~=ic);
    nc = length(ic2);
    du = zeros(nc,4);
    dX = zeros(nc,2);
    for j=1:nc
        du(j,:) = phi(ic2(j),:) - phi(ic,:);
        disp([num2str(ic),',',num2str(ic2(j)),': du(',num2str(j),...
             ',:) = ',num2str(du(j,:))]);
        %dX(j,:) = XYC(ic2(j),:) - XYC(ic,:);
    end
    %dphidX(:,:,ic) = (dX\du)';
    dphidX(:,:,ic) = (dXinv(:,1:nc,ic)*du)';
end
toc