function [dH ] = RHS(t , fseta ,x, ph  , h , mi ) 
gi = 9.81; % accelaration of gravity

Nx = length(x); % the number of grid points, including boundaries 
eta = fseta(1:Nx,2);
psi = fseta(1:Nx,1); % free surface potential
dx = x(2)-x(1);
psix = gradientp(psi,dx);
etax = gradientp(eta,dx);
% 
deta = -psix(1:Nx).*etax(1:Nx)+(etax(1:Nx).^2+1).*(ph(1:Nx)/h + mi*psi(1:Nx))  ;
dpsi = -gi*eta(1:Nx) - 0.5*psix(1:Nx).^2 + 0.5*(etax(1:Nx).^2+1).*(ph(1:Nx)/h + mi*psi(1:Nx)).^2 ;

dH = [dpsi  deta];



