% -------------------------------------------------------------------------
%
%                         FV2D Post Processor                           
%           Plots data generated & saved from FV2D Euler solver
%                       
%                     AA 214C, Stanford University
%                       Jacob Crabill  4/21/2014
% -------------------------------------------------------------------------
clear all;
close all;

addpath('./data/');

%% ------------------------------------------------------------------------
% Data Import
% -------------------------------------------------------------------------
caseName = 'ViscBump_Fine';
caseIter = 10000;

dataFile = [caseName,'_',num2str(caseIter)];
load(dataFile);

% Plotting Limits
% xmin = 0;
% xmax = 10;
% ymin = 0;
% ymax = 2.5;

%% ------------------------------------------------------------------------
% Variable Setup
% -------------------------------------------------------------------------
ic = (1:n_cells)';
x = reshape(xv(c2v(ic,:),1),n_cells,max_nf)';
y = reshape(xv(c2v(ic,:),2),n_cells,max_nf)';
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

rho = W(ic,1);
u = W(ic,2)./rho;
v = W(ic,3)./rho;
E = W(ic,4);
p = (gam-1)*(E - 0.5*rho.*(u.^2+v.^2));
asq = gam*p./rho;
T = asq/(gam*R)^2;
M = sqrt((u.^2+v.^2)./asq);

%% ------------------------------------------------------------------------
% Solution Plots
% -------------------------------------------------------------------------
% % Density
% h1 = figure(1); grid on,xlabel('x'),ylabel('y'),
% xlim([xmin,xmax]),ylim([ymin,ymax]);
% patch(x,y,p','EdgeColor','none');
% title('Density'), colorbar, drawnow;
% print(h1,'-dpng',['./plots/',dataFile,'_rho']);
% 
% % Pressure
% h2 = figure(2); grid on,xlabel('x'),ylabel('y'),
% xlim([xmin,xmax]),ylim([ymin,ymax]);
% patch(x,y,p','EdgeColor','none');
% title('Pressure'), colorbar, drawnow;
% print(h2,'-dpng',['./plots/',dataFile,'_P']);
% 
% % % Mach Number
% h3 = figure(3); grid on,xlabel('x'),ylabel('y'),
% xlim([xmin,xmax]),ylim([ymin,ymax]);
% patch(x,y,M','EdgeColor','none');
% title('Mach Number'), colorbar, drawnow;
% % print(h3,'-dpng',['./plots/',dataFile,'_M']);
% % 
% % % X Velocity
% h4 = figure(4); grid on,xlabel('x'),ylabel('y'),
% xlim([xmin,xmax]),ylim([ymin,ymax]);
% patch(x,y,u','EdgeColor','none');
% title('X-Velocity'), colorbar, drawnow;
% % print(h4,'-dpng',['./plots/',dataFile,'_U']);
% % 
% % % Y Velocity
% h5 = figure(5); grid on,xlabel('x'),ylabel('y'),
% xlim([xmin,xmax]),ylim([ymin,ymax]);
% patch(x,y,v','EdgeColor','none');
% title('Y-Velocity'), colorbar, drawnow;
% % print(h5,'-dpng',['./plots/',dataFile,'_V']);

%% ------------------------------------------------------------------------
% Residual
% -------------------------------------------------------------------------
h6 = figure(6);
semilogy(res,'linewidth',2); xlabel('Iteration'); ylabel('residual');
grid on; xlim([0 iter]); ylim([1E-6 1]);
legend('\rho','\rhou','\rhov','E');
title(['Residual, iteration ',num2str(iter)]); drawnow;
print(h6,'-dpng',['./plots/',dataFile,'_res']);

%% ------------------------------------------------------------------------
% Cp over bottom wall
% -------------------------------------------------------------------------
rho_inf = p_inf/(R*T_inf);         % kg/m^3    | Density
a_inf = sqrt(gam*R*T_inf);         % m/s       | Speed of sound
n_inf = [cos(theta);sin(theta)];   % --        | Flow direction
u_inf = M_inf*a_inf*n_inf(1);      % m/s       | x-direction velocity
v_inf = M_inf*a_inf*n_inf(2);      % m/s       | y-direction velocity
E_inf = p_inf/(gam-1) + rho_inf*(M_inf*a_inf)^2/2; % J/m^3 = kg/m-s^2 | Energy
Tt_inf = T_inf*(1+(gam-1)/2*M_inf^2);  % K   | Total temperature
pt_inf = p_inf*(1+(gam-1)/2*M_inf^2)^(gam/(gam-1));  % Pa  | Total pressure

ie = [ie_slip_wall;ie_adiabatic_wall];
ie = ie(XYE(ie,2)<.1);
ie = ie(XYE(ie,1) > -.2);
ie = ie(XYE(ie,1) < .3);
ic = e2c(ie,1);
ic2 = e2c(ie,2);

un = phi(ic,2).*unorm(ie,1)+phi(ic,3).*unorm(ie,2);
c = gam*phi(ic,4)./phi(ic,1);

dX = 2*(XYE(ie,:)-XYC(ic,:));
pwall = phi(ic,4) + squeeze(dphidX(4,1,ic)).*dX(:,1) + squeeze(dphidX(4,2,ic)).*dX(:,2);
Cp = (pwall - p_inf)/(0.5*rho_inf*u_inf^2);

XCP = sortrows([XYE(ie,1),-Cp]);
X = XCP(:,1);
nCp = XCP(:,2);

h7 = figure(7);
plot(X,nCp,'-k*'); grid on; xlabel('x/c'); ylabel('-Cp');
title('Neg. Pressure Coeff. - Mach 0.908, 6% thickness circular arc');
print(h7,'-dpng',['./plots/',dataFile,'_Cp']);

%% ------------------------------------------------------------------------
% Mach Contours Near Bump
% -------------------------------------------------------------------------
h8 = figure(8);
x1 = linspace(-.1,.2,250);
y1 = linspace(0,.15,250);
[X,Y]=meshgrid(x1,y1);
ic = 1:n_cells;
M = sqrt(sum(phi(ic,[2,3]).^2,2)./(1.4*phi(ic,4)./phi(ic,1)));
F = TriScatteredInterp(XYC(1:n_cells,1),XYC(1:n_cells,2),M);
M2 = F(X,Y);

xc = .05;
yc = -.2053;
r = .2113;

contourf(X,Y,M2);
hold on;
rectangle('Position',[xc-r,yc-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');         
xlabel('X (m)'); ylabel('Y (m)'); colorbar;
title('Mach contours in region of bump');
print(h8,'-dpng',['./plots/',dataFile,'_M_contour']);

%% ------------------------------------------------------------------------
% Velocity Contours Near Bump
% -------------------------------------------------------------------------
% Velocity Magnitude
h11 = figure(11);
x1 = linspace(.08,.15,100);
y1 = linspace(0,.03,100);
[X,Y]=meshgrid(x1,y1);
ic = 1:n_cells;
magU = sqrt(sum(phi(ic,[2,3]).^2,2));
F = TriScatteredInterp(XYC(1:n_cells,1),XYC(1:n_cells,2),magU);
magU = F(X,Y);

xc = .05;
yc = -.2053;
r = .2113;

contourf(X,Y,magU);
hold on;
rectangle('Position',[xc-r,yc-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');         
xlabel('X (m)'); ylabel('Y (m)'); colorbar;
title('Velocity Magnitude Contours near Trailing Edge');
print(h11,'-dpng',['./plots/',dataFile,'_magU_contour']);

% X-Velocity
x1 = linspace(.08,.12,100);
y1 = linspace(0,.015,100);
[X,Y]=meshgrid(x1,y1);
ic = 1:n_cells;

h12 = figure(12);
F = TriScatteredInterp(XYC(ic,1),XYC(ic,2),phi(ic,2));
U = F(X,Y);

xc = .05;
yc = -.2053;
r = .2113;

contourf(X,Y,U);
hold on;
rectangle('Position',[xc-r,yc-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');         
xlabel('X (m)'); ylabel('Y (m)'); colorbar;
title('X-Velocity Contours near Trailing Edge');
print(h12,'-dpng',['./plots/',dataFile,'_U_contour']);

%% ------------------------------------------------------------------------
%  Velocity Vectors
% -------------------------------------------------------------------------
h9 = figure(9);
Y = XYC(:,2);
X = XYC(:,1);
ic = find(Y<=.015);
ic = ic(X(ic)>.08);
ic = ic(X(ic)<.12);
x = XYC(ic,1);
y = XYC(ic,2);
u = phi(ic,2);
v = phi(ic,3);

L = .1;
deltay = .06*L;
xc = L/2;
yc = deltay/2 - L^2/(8*deltay);
r = deltay - yc;

xbump = XYE(ie_adiabatic_wall,1);
ybump = XYE(ie_adiabatic_wall,2);
XY = sortrows([xbump,ybump]);
xbump = XY(:,1);
ybump = XY(:,2);

%rectangle('Position',[xc-r,yc-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');
plot(xbump,ybump,'k-');
xlim([.08,.12]); ylim([0,.015]);
hold on; grid on;
quiver(x,y,u,v,.3,'k','filled','LineWidth',2,'MaxHeadSize',.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Velocity vectors near trailing edge');
print(h9,'-dpng',['./plots/',dataFile,'_TE_Vectors']);

% Zoom near TE
h10 = figure(10);
ic = find(Y<=.01);
ic = ic(X(ic)>.08);
ic = ic(X(ic)<.1);
x = XYC(ic,1);
y = XYC(ic,2);
u = phi(ic,2);
v = phi(ic,3);

plot(xbump,ybump,'k-');
%rectangle('Position',[xc-r,yc-r,2*r,2*r],'Curvature',[1,1],'FaceColor','w');
xlim([.08,.1]); ylim([0,.01]);
hold on; grid on;
quiver(x,y,u,v,.3,'k','filled','LineWidth',2,'MaxHeadSize',.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Velocity vectors near trailing edge');
print(h10,'-dpng',['./plots/',dataFile,'_TE_Vectors_Zoom']);