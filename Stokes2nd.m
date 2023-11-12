function [ eta, psi ] = Stokes2nd( a,k, h, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gi = 9.81;
om = sqrt(gi*k*tanh(k*h));
mu2 = 0.5*k*coth(k*h)*(1+3/(2*sinh(k*h)^2));
nu1 = om/k/sinh(k*h);
nu2 = (3/8)*om/sinh(k*h)^4;
eta = a*cos(k*x) + mu2*a*a*cos(2*k*x);
psi = nu1*a*cosh(k*(eta+h)).*sin(k*x) + nu2*a*a*cosh(2*k*(eta+h)).*sin(2*k*x);
eta = eta';
psi = psi';
end

