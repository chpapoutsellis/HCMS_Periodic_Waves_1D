clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCMS, periodic boundary conditions, 4th order Finite-Difference, RK4    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gi = 9.81; % accelaration of gravity
h0  = 1;    % depth [m]

%----------------------- INITIAL CONDITIONS ------------------------------%
a = 0.1; % amplitude
L = 8;   % wave length

x = linspace(0,L,128);
Nx = length(x);

h = ones(length(x),1); % depth

dx = x(2)-x(1);
k0 = 2*pi/L; % wave number for the Stokes wave
[ eta0, fi0 ] = Stokes2nd( a,1*k0, h0, x );

% perturbed Stokes wave
delta = 0;
epsilon = 0;
eta0p = eta0 + epsilon*delta.*eta0;
fi0p = fi0 + epsilon*delta.*fi0;


% Initialize figure
close all
figure('color','w','position', [287         522        1319         420])
subplot(2,1,1)
hpert_eta = line(x,eta0p);
ylim([-1.5*a 1.5*a])
ylabel('$\eta$','interpreter','latex','fontsize',20)
grid on
box on
xlim([x(1) x(end)])
subplot(2,1,2)
hpert_psi = plot(x,fi0p,'r');hold on;
xlim([x(1) x(end)])
xlabel('$x$','interpreter','latex','fontsize',20)
ylabel('$\psi$','interpreter','latex','fontsize',20)
hold on
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CMS params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nm = 5+1;  % number of modes in the HCMS %
om = sqrt(gi*k0*tanh(k0*h0)); % finite-depth linear frequency corresponding to k
mi0 = om*om/gi;               % Choice of the parameter \mu_0 in the vertical functions
% h0 is the depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = 2*pi/om; % period from linear dispersion
c_sh = sqrt(gi*h0); 
dt = 0.75*dx/c_sh;
ti = 0:dt:10*T; % Time grid
dt = ti(2)-ti(1);
NT = length(ti); 
COUR = c_sh*dt/dx;

fi = zeros(Nx,Nm,NT);
eta = zeros(Nx,NT);
fs = eta;
Ham = zeros(1,NT);
Mass = Ham;
Qx = Ham;

%%%%%%%%%%%%%%%%%%%%% RK4$  time integrator  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fi0 = real(fi0);
U0 = [fi0p eta0p];


% Initial fields and conserved quantities
etx = gradientp(eta0,dx);
fix = gradientp(fi0,dx);
[ph] = dtncp_flat5_sq(0, U0 , Nm-1 , x, h, mi0, h0);
Hamd0 = 0.5*(fi0.*(-etx.*fix + (etx.^2 + 1).*(ph/h0 + mi0*fi0))+ gi*eta0.^2);
Ham0 = trapz(x,Hamd0);  % Energy
M0 = trapz(x(1:Nx),eta0(1:Nx)); % Mass
Qx0 = trapz(x,eta0.*gradientp(fi0,dx)); % Horizontal impulse
eta(:,1) = U0(1:Nx,2); % free-surface elevation
fs(:,1) = U0(1:Nx,1);  % free-surface potential
    

% print timestep and errors in Hamiltonian, Mass and Momentum
fprintf('   it|   Err{H}      |  Err{M}      |  Err{Q}      \n');

for it = 2:length(ti)
    tic;
    etx = gradientp(U0(Nx+1:2*Nx),dx);
    fix = gradientp(U0(1:Nx),dx);
    
    % Solution of the CMS
    [ph] = dtncp_flat5_sq(0, U0 , Nm-1 , x,h, mi0, h0);

    % Calculation of conserved quantities
    Hamd = 0.5*(U0(1:Nx,1).*(-etx.*fix + (etx.^2 + 1).*(ph/h0 + mi0*U0(1:Nx,1)))+ gi*U0(1:Nx,2).^2);
    Ham(it) = trapz(x,Hamd);
    Mass(it) = trapz(x,U0(1:Nx,2));
    Qx(it) = trapz(x,U0(1:Nx,2).*fix);

    % Right hand side of the HCMS
    [ f1 ] = RHS(ti(it), U0 , x, ph  , h0 , mi0);

    
    [ph] = dtncp_flat5_sq(0, U0 + 0.5*dt*f1 , Nm-1 , x,h, mi0, h0);
    [ f2 ] = RHS(ti(it)+dt/2, U0 + 0.5*dt*f1, x ,ph , h0 , mi0);

    [ph] = dtncp_flat5_sq(0, U0 + 0.5*dt*f2 , Nm-1 , x,h, mi0, h0);
    [ f3 ] = RHS(ti(it)+dt/2, U0 + 0.5*dt*f2, x ,ph  , h0 , mi0 );

    [ph] = dtncp_flat5_sq(0, U0 + dt*f3,  Nm-1 , x,h, mi0, h0);
    [ f4 ] = RHS(ti(it)+dt, U0 + dt*f3, x ,ph  , h0 , mi0);

    Ut = U0 + (1/6)*dt*( f1 + 2*f2 + 2*f3 + f4);
    U0 = Ut;
    t1 = toc;

    
    eta(:,it) = U0(1:Nx,2);
    fs(:,it) = U0(1:Nx,1);
    CompTime = (it)*t1;
    RealTime = ti(it);

    fprintf('    %d|   %1.5f      |  %1.5f      |  %1.5f      \n', it ,(abs(Ham(it)-Ham0)/Ham0), (abs(Mass(it))),(abs(Qx(it)-Qx0)/Qx0));
    if mod(it,20)==0
        subplot(2,1,1)
        set(hpert_eta,'Xdata',x,'Ydata',U0(1:Nx,2))
        subplot(2,1,2)
        set(hpert_psi,'Xdata',x,'Ydata',U0(1:Nx,1))
        drawnow
    else
    end
end

%%%%%%%%%%%%%%%%%% End of RK4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Space-time plot
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
fs= 16; % fontsize
close all
[XX,TT] = meshgrid(x/L,ti/T);
figure('color','w','position',[519   557   735   293])
pcolor(XX,TT,real(eta)')
clim([-a a])
shading interp
colormap(parula)
colorbar
xlabel('$x/L$','Interpreter','latex')
ylabel('$t/T$','Interpreter','latex')
set(gca,'fontsize',fs)



