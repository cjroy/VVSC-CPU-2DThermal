% clear all
close all
clc

%============Inputs============

%**System Paremeters**
TDP = 100.0;           % Thermal Design Power, W
NumFins = 5;           % Number of heat sink fins
w = 8.0/100;           % Total length of domain parallel to flow direction, m
L = 10.0/100;           % Total length of domain perpendicular flow direction (including L_CPU), m
L_CPU = 2.0/100;       % Length of CPU
thick = 0.4/100;       % Thickness of heat sink fin, m

%**Material Parameters**
k_Metal = 160.;               % Thermal conductivity of Metal, W/(m*K)
density_Metal = 2700.;        % Density of metal, kg/m^3
specific_heat_Metal = 895.;   % Specific heat of metal, J/(kg*K)
alpha_Metal = k_Metal/(density_Metal*specific_heat_Metal);   % Thermal diffusivity of Metal, m^2/s

k_CPU = 160.00;               % Thermal conductivity of CPU, W/(m*K)
density_CPU = 2700.;          % Density of CPU, kg/m^3
specific_heat_CPU = 895.;     % Specific heat of CPU, J/(kg*K)
alpha_CPU = k_CPU/(density_CPU*specific_heat_CPU); % Thermal diffusivity of CPU, m^2/s

%**Fluid Properties**
Vel_air = 10.0;            % m/s (Set Vel_air=0 to use hbar, otherwise Vel_air is used to determine convective cooling coefficient)
hbar = 15.0;              % Averaged convective cooling coeff., W/(m^2*K) (Not used unless Vel_air=0)
Tinf = 300.0;             % Averaged air temperature, K

%**Discretization Parameters**
imax = 17;                % Number of nodes in x (normal to flow direction) 
jmax = 17;                % Number of nodes in x (tangential to flow direction) 
itermax = 2e7;            % Maximum allowable number of iterations
convtol = 1.e-8;         % Iterative convergence tolerance (relative to fifth iteration)

%**Solution Variable**
T = 350.*ones(imax,jmax); % Temperature of heat sink (initialized to 350 K), K


%%**Calculate Solutin**
tic
% [T, x, y, L2conv, history] = heatcondsolve2d(T,TDP,NumFins,Lx,Ly,thick,k_Metal,alpha_Metal,hbar,Tbar,imax,jmax,itermax,convtol);
[T, x, y, L2conv, history] =heatcondsolve2d(T,TDP,NumFins,w,L,thick,k_Metal,alpha_Metal,k_CPU,alpha_CPU,L_CPU,Vel_air,hbar,Tinf,imax,jmax,itermax,convtol);

toc

MaxTemp = max(max(T)) - 273.15;
fprintf('max Temp: %23.18f\n',MaxTemp)

%calculating integral of temperature


% NOTE: you may want to suppress figure output for nondeterministic simulations
figure(1)
subplot(2,1,1)
semilogy(history(:,1),history(:,2))
title('Iterative convergence')

subplot(2,1,2)
surface(y,x,T-273.15)
view([0,0,90])
title('Temperature (deg C)')
xlabel('W')
ylabel('L')


figure(2) % NOTE: Figure2 plots the temperature profiles at center of plate
subplot(2,1,1)
plot( x(:,(jmax-1)/2+1),T(:,(jmax-1)/2+1)-273.15 )
title('Temperature (deg C) vs. x')

subplot(2,1,2)
plot( y((imax-1)/2+1,:),T((imax-1)/2+1,:)-273.15 )
title('Temperature (deg C) vs. y')
