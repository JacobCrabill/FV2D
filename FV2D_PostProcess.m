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
caseName = 'FlatPlate_MC';
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
% % Mach Number
h3 = figure(3); grid on,xlabel('x'),ylabel('y'),
xlim([xmin,xmax]),ylim([ymin,ymax]);
patch(x,y,M','EdgeColor','none');
title('Mach Number'), colorbar, drawnow;
% print(h3,'-dpng',['./plots/',dataFile,'_M']);
% 
% % X Velocity
h4 = figure(4); grid on,xlabel('x'),ylabel('y'),
xlim([xmin,xmax]),ylim([ymin,ymax]);
patch(x,y,u','EdgeColor','none');
title('X-Velocity'), colorbar, drawnow;
% print(h4,'-dpng',['./plots/',dataFile,'_U']);
% 
% % Y Velocity
h5 = figure(5); grid on,xlabel('x'),ylabel('y'),
xlim([xmin,xmax]),ylim([ymin,ymax]);
patch(x,y,v','EdgeColor','none');
title('Y-Velocity'), colorbar, drawnow;
% print(h5,'-dpng',['./plots/',dataFile,'_V']);

%% ------------------------------------------------------------------------
% Residual
% -------------------------------------------------------------------------
% h6 = figure(6);
% semilogy(res,'linewidth',2); xlabel('Iteration'); ylabel('residual');
% grid on; xlim([0 iter]); ylim([1E-6 1]);
% legend('\rho','\rhou','\rhov','E');
% title(['Residual, iteration ',num2str(iter)]); drawnow;
% print(h6,'-dpng',['./plots/',dataFile,'_res']);

%% ------------------------------------------------------------------------
% Cp over bottom wall
% -------------------------------------------------------------------------
% rho_inf = p_inf/(R*T_inf);         % kg/m^3    | Density
% a_inf = sqrt(gam*R*T_inf);         % m/s       | Speed of sound
% n_inf = [cos(theta);sin(theta)];   % --        | Flow direction
% u_inf = M_inf*a_inf*n_inf(1);      % m/s       | x-direction velocity
% v_inf = M_inf*a_inf*n_inf(2);      % m/s       | y-direction velocity
% E_inf = p_inf/(gam-1) + rho_inf*(M_inf*a_inf)^2/2; % J/m^3 = kg/m-s^2 | Energy
% Tt_inf = T_inf*(1+(gam-1)/2*M_inf^2);  % K   | Total temperature
% pt_inf = p_inf*(1+(gam-1)/2*M_inf^2)^(gam/(gam-1));  % Pa  | Total pressure
% 
% ie = ie_slip_wall;
% ie = ie(XYE(ie,2)<1);
% ie = ie(XYE(ie,1) > 0);
% ie = ie(XYE(ie,1) < 1);
% ic = e2c(ie,1);
% ic2 = e2c(ie,2);
% 
% un = phi(ic,2).*unorm(ie,1)+phi(ic,3).*unorm(ie,2);
% c = gam*phi(ic,4)./phi(ic,1);
% 
% dX = 2*(XYE(ie,:)-XYC(ic,:));
% pwall = phi(ic,4) + squeeze(dphidX(4,1,ic)).*dX(:,1) + squeeze(dphidX(4,2,ic)).*dX(:,2);
% Cp = (pwall - p_inf)/(0.5*rho_inf*u_inf^2);
% 
% h7 = figure(7);
% plot(XYE(ie,1),-Cp,'-k*'); grid on; xlabel('x/c'); ylabel('-Cp');
% title('Neg. Pressure Coeff. - Mach 0.8, 6% thickness bump');
% print(h7,'-dpng',['./plots/',dataFile,'_Cp']);

%% ------------------------------------------------------------------------
% Velocity Profile at x = L/2
% -------------------------------------------------------------------------
h8 = figure(8);
x = .0005;
ic1 = find(XYC(:,1)==x);
y2 = .0000775;
ic = ic1(XYC(ic1,2)<=y2);
ic = ic(XYC(ic,2)>=0);
u = phi(ic,2);

% -----------------------------
% Blasius Solution
% -----------------------------
% Solving the ODE with the refined initial condition on F_alphaalpha
bc = .46975;
Ue = max(u);
ye = XYC(ic(u==Ue),2);
[alpha,f]=ode45(@Blasius,[0 6],[0 0 bc]);

% Freestream Conditions
mu_inf = 1.8865E-5;
rho_inf = p_inf/(R*T_inf);
nu = mu_inf / rho_inf;

% Map Blasius solution to dimensional solution
y = alpha/sqrt(Ue/(2*nu*x));
UUe = f(:,2);

plot(Ue*UUe,y,'r--','LineWidth',2);
hold on;
plot(u,XYC(ic,2),'k-*','LineWidth',2);
legend('Blasius','Computed');
title('Boundary Layer Velocity Profile at X = L/2');
grid on; xlabel('x-velocity (m/s)'); ylabel('y position (m)');
%print(h8,'-dpng',['./plots/',dataFile,'_U_Profile']);

%% ------------------------------------------------------------------------
%  Velocity Vectors
% -------------------------------------------------------------------------
h9 = figure(9);
Y = XYC(:,2);
ic = find(Y<2E-4);
x = XYC(ic,1);
y = XYC(ic,2);
u = phi(ic,2);
v = phi(ic,3);
quiver(x,y,u,v,.5,'k','filled','LineWidth',2,'MaxHeadSize',.002);
grid on;

[X,Y]=meshgrid(2E-5:.001/200:.001,0:.001/500:7E-5);
F = TriScatteredInterp(XYC(:,1),XYC(:,2),phi(:,2));
U = F(X,Y);

mu_inf = 1.8865E-5;
rho_inf = p_inf/(R*T_inf);
Uinf = 106.375;
hold on;
contour(X,Y,U,.99*Uinf,'r','LineWidth',2)

% % Plot Blasius profile as well
% x2 = linspace(0,.001,100);
% delta = 5*x2./sqrt(rho_inf*Uinf*x2/mu_inf);%1.9949E-3*sqrt(x2);
% plot(x2,delta,'b--','LineWidth',2);
% legend('Vector Field','Computed 95% line','Blasius 99% line');

xlim([xmin, xmax]); ylim([0, 2E-4]);
xlabel('X (m)'); ylabel('Y (m)');
title('Velocity Vectors near Boundary Layer, and .95U_e Contour');
print(h9,'-dpng',['./plots/',dataFile,'_U_Vectors']);

%% ------------------------------------------------------------------------
%  Friction Coefficient
% -------------------------------------------------------------------------
rho_inf = p_inf/(R*T_inf);
mu_inf = 1.8865E-5;
Uinf = 106.375;

pwall = phiL(ie_adiabatic_wall,4);
rhowall = phiL(ie_adiabatic_wall,1);
Twall = pwall./(R*rhowall);

T_s = 110.56;
T_ref = 273.11;
mu_ref = 1.716E-5;
mu = mu_ref*(Twall./T_ref).^(3/2).*((T_ref+T_s)./(Twall+T_s));

dudy = tauL(ie_adiabatic_wall,2);
dvdx = tauL(ie_adiabatic_wall,3);

x = XYE(ie_adiabatic_wall,1);

tauwall = mu.*(dudy+dvdx);
Cf1 = tauwall/(1/2*rho_inf*Uinf^2);
Cf2 = .664./sqrt(rho_inf*Uinf*x/mu_inf);
plot(x,Cf1,'k-*');
grid on; hold on;
plot(x,Cf2,'r-*');
xlabel('X (m)'); ylabel('C_f');
legend('Computed Result','Blasius Solution');
title('Friction Coefficien on Wall');

%% Misc 1
[X,Y] = meshgrid(0:.001/200:.001);
M = sqrt(sum(phi(:,[2,3]).^2,2)./(1.4*phi(:,4)./phi(:,1)));
F = TriScatteredInterp(XYC(:,1),XYC(:,2),M);
m = F(X,Y);
contourf(X,Y,m,10)
xlim([0,.001]);ylim([0,2E-4]);
xlabel('X (m)'), ylabel('Y (m)')
title('Mach contours near boundary layer');

%% Misc 2
[X,Y] = meshgrid(0:.001/200:.001);
U = phi(:,3);
F = TriScatteredInterp(XYC(:,1),XYC(:,2),U);
u = F(X,Y);
contourf(X,Y,u,15); colorbar
xlim([0,.001]);ylim([0,2E-4]);
xlabel('X (m)'), ylabel('Y (m)')
title('Y-velocity contours near boundary layer (m/s)');