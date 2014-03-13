function [T, x, y, L2conv, history] = heatcondsolve2d(T,TDP,NumFins,w,L,thick,k_Metal,alpha_Metal,k_CPU,alpha_CPU,L_CPU,Vel_air,hbar,Tinf,imax,jmax,itermax,convtol)
% Solves the 2D heat equation with a specified energy flux to approximate 
% the base temperature of a heat sink with fins in a laminar fluid medium.
%
%
%                  ___________________________
%                 /                          /|
%                ---------------------------/ |
%                |                          | |
%                |                          | |   x (L)
%   -----> V     |          FIN             | |   ^
%                |                          | |   |  z (thick)
%                |                          | |   | /
%                |                          | |   |/
%                ----------------------------/|    -----> y (w)
%                |          CPU             | |
%                ----------------------------/
%
%
% [Tbase, T, x, L2conv, history] =
% heatcondsolve(T,TDP,NumFins,L,thick,L_CPU,k_Metal,alpha_Metal,k_CPU,alpha_CPU,hbar,Tinf,imax,itermax,convtol)
% 
%
% Inputs:
%                     T : Initial Temperature (matrix [1 x imax] or [imax x 1] )
%                   TDP : Thermal Design Power, W 
%               NumFins : Number of heat sink fins 
%                     w : Length of fins in flow direction, m 
%                     L : Length of fins perpendicular to flow direction, m 
%                 thick : Thickness of heat sink fin, m 
%               k_Metal : Thermal conductivity of Metal, W/(m*K)
%           alpha_Metal : Thermal diffusivity of Metal, m^2/s
%                 k_CPU : Thermal conductivity of CPU, W/(m*K)
%             alpha_CPU : Thermal diffusivity of CPU, m^2/s
%                 L_CPU : Length of CPU
%               Vel_air : Air velocity to determine convective cooling coefficient, W/(m^2*K)
%                  hbar : Average convective cooling coefficient (Set Vel_air=0 to use hbar), W/(m^2*K)
%                  Tinf : Air temperature, K
%                  imax : number of nodes in x
%                  jmax : number of nodes in y
%               itermax : maximum number of iterations
%               convtol : Normalized convergence tolerance
% 
% Outputs: 
%                     T : Final temperature
%                     x : x node locations
%                     y : y node locations
%                L2conv : Final normalized residual
%               history : Iterative residual history  of every 100th iteration [iter #, residual]


%% Initialize variables
history=[0,0];          % Initialize history vector
L2normalize=1;          % Normalizing value for iterative residuals
insul = zeros(size(T)); % Insulation blanking array
k=insul;                % Thermal conductivity
alpha=insul;            % Thermal diffusivity
RHS=zeros(size(T));     % Initialize RHS vector used in iterative scheme
history_output=1000;    % How frequently to store iterative residuals

%% Error checking
if any(size(T)~= [imax,jmax])
    error('Temperature vector must match specified number of nodes!')
end


%% Set convective coefficient
mu_air=1.95e-5; %Ns/m^2
Cp_air=1008;   %J/kgK
K_air=2.8e-2;   %W/mK
rho_air=1.05;   %kg/m^3
Pr=mu_air*Cp_air/K_air; %Prandtl number
hfun = @(x,y)0.332*K_air*Pr^(1/3)*(rho_air*Vel_air./(mu_air*(y+0.5*thick))).^(1/2);           % Convective coefficient as a function of y, W/(m^2*K)




%% Set geometry
x=linspace(0,L,imax);
y=linspace(0,w,jmax);
[x,y]=meshgrid(x,y);
x=x';y=y';
dx = L/(imax-1);
dy = w/(jmax-1);

%% Set thermal conductivity and thermal diffusivity arrays
alpha=alpha_Metal*ones(imax,jmax);
k=k_Metal*ones(imax,jmax);
if Vel_air>0
    h=hfun(x,y);
else
    h=hbar*ones(size(x));
end


if L_CPU~=0
    for i = 1:imax
        k(i,:) = k_CPU;
        alpha(i,:) = alpha_CPU;
        insul(i,:) = 0.;
        if x(i,1) > L_CPU 
            insul(i,:) = 1.0;
            k(i,:) = k_Metal;
            alpha(i,:) = alpha_Metal;
        end 
        if abs(x(i,1) - L_CPU)<1.e-6
            insul(i,:) = 0.5;
            k(i,:) = 0.5*(k_Metal + k_CPU);
            alpha(i,:) = 0.5*(alpha_Metal + alpha_CPU);
        end
    end
else
    insul=ones(size(T));
end
    


%% Set time step array (local time stepping)
dt = 0.225*dx*dy./alpha;


%% Iterative solution loop
i=2:imax-1;
j=2:jmax-1;
src_constant = 2*alpha.*h./(k.*thick).*insul;
hist_cnt = 1;
for iter = 1:itermax

    % Compute Right Hand Side
    d2Tdx2 = ( T(i+1,j) - 2*T(i,j) + T(i-1,j) )/(dx^2);
    d2Tdy2 = ( T(i,j+1) - 2*T(i,j) + T(i,j-1) )/(dy^2);
    src = src_constant(i,j).*(T(i,j)-Tinf);
    RHS(i,j) = alpha(i,j).*(d2Tdx2+d2Tdy2) - src;
    L2conv = sqrt( sum(sum(RHS(i,j).^2))/((imax-2)*(jmax-2)) );

    % Update Temperatures in Interior
    T(i,j) = T(i,j) + dt(i,j).*RHS(i,j);

    %**Update Boundary Conditions**
    %@y=0, j=1
    ibc=2:imax-1;
    jbc=1;
    T(ibc,jbc)=(4*T(ibc,jbc+1)-T(ibc,jbc+2))./3;

    %@y=w, j=jmax
    ibc=2:imax-1;
    jbc=jmax;
    T(ibc,jbc)=(4*T(ibc,jbc-1)-T(ibc,jbc-2))./3;

    %@x=L, i=imax
    ibc=imax;
    jbc=1:jmax;
    T(ibc,jbc)=(4*T(ibc-1,jbc)-T(ibc-2,jbc))./3;

    %@x=0, i=1
    ibc=1;
    jbc=1:jmax;
    psi2=2*TDP*dx./(NumFins*k(ibc,jbc)*thick*w);
    T(ibc,jbc)=(4*T(ibc+1,jbc)-T(ibc+2,jbc)+psi2)./3;

    % Normalize residuals
    if (iter==5); L2normalize=L2conv; end
    L2conv=L2conv/L2normalize;
    
    % Store iterative residuals every 100 iterations
    if mod(iter,history_output)==0;
        history(hist_cnt,:) = [iter,L2conv];
        hist_cnt=hist_cnt+1; 
%         fprintf('%7.0f %12.6e %23.15e\n',iter,L2conv,max(max(T))); %To monitor convergence. Slows the convergence time because of printing
    end

    % Check for iterative convergence
    if (L2conv < convtol) && (iter>2); break ; end
        
end

%% Convergence checking
if (L2conv > convtol)
    fprintf('Warning! Solution not Converged. Imax = %4.0f, Residual = %5.3e\n',imax, L2conv)
end


end
