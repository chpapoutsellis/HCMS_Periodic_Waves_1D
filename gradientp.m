function f = gradientp(y,dx)
% Calculates the finite difference approximation of the derivative of the
% periodic function y defined on x.
 

N = length(y);
f = zeros(N,1);

for i=3:N-2
f(i) = (-y(i+2)+8*y(i+1)-8*y(i-1)+y(i-2))/12/dx;
end
f(1) = ( -y(3) + 8*y(2) - 8*y(end-1) + y(end-2) )/12/dx;
f(2) = ( -y(4) + 8*y(3) - 8*y(1) + y(end-1) )/12/dx;
f(N) = (-y(3) + 8*y(2) - 8*y(N-1) + y(N-2))/12/dx;
f(N-1) = (-y(2) + 8*y(N) - 8*y(N-2) + y(N-3))/12/dx;















