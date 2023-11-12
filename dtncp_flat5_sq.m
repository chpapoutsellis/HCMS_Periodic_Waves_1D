function [ph]=dtncp_flat5_sq(t, U , Nm , x, h, mi0, h0 ) 
% Solution of the CMS with a 4th-order finite-difference method and
% periodic boundary conditions
% INPUTS
%   Nm      : Number of modes
% As, Bs, Cs: matrix coefficients
%   x       : horizontal grid
%   U       : (psi;eta)
% mi0,h0    : parameters for the vertical functions
%   h       : depth

% OUTPUT:
%   ph      : ph_{-2} modal function


Nx = length(x); % the number of grid points, including boundaries 
psi = U(1:length(x));
dx = x(2)-x(1);

% Calculation of the matrix coefficients of the CMS
[ As, Bs, Cs  ] = coef_ABC_flat(  U(1:Nx,2) , h , Nm+1 , mi0 ,dx ,h0 );

% Construction and solution of the linear system
Am = zeros((Nm-1)*Nx,Nm*Nx);
Bv = zeros(Nm*Nx,1);
   
for m=1:Nm-1
    A1 = As(:,:,1);
    B1 = Bs(:,:,1);
    C1 = Cs(:,:,1);
    AA = A1;
    BB = B1;
    CC = C1;

    for n=1:Nm
        Am((m-1)*Nx+1,(n-1)*Nx+1) = -30*AA(m,n)/dx/(dx*12) + CC(m,n);
        Am((m-1)*Nx+1,(n-1)*Nx+2) = 16*AA(m,n)/dx/(dx*12) + 8*(BB(m,n))/12/dx ;
        Am((m-1)*Nx+1,(n-1)*Nx+3) = (-(AA(m,n))/dx/(dx*12)) - ((BB(m,n))/12/dx) ;
        Am((m-1)*Nx+1,(n-1)*Nx+Nx-1) = ((16*AA(m,n))/dx/(dx*12)) - (8*(BB(m,n))/12/dx);
        Am((m-1)*Nx+1,(n-1)*Nx+Nx-2) = (-(AA(m,n))/dx/(dx*12)) + ((BB(m,n))/12/dx);
    end
    A1 = As(:,:,end);
    B1 = Bs(:,:,end);
    C1 = Cs(:,:,end);
    AA = A1;
    BB = B1;
    CC = C1;
    for n=1:Nm
        Am((m)*Nx,(n-1)*Nx+2) = ((16*AA(m,n))/dx/(dx*12)) + (8*(BB(m,n))/12/dx);
        Am((m)*Nx,(n-1)*Nx+3) = (-(AA(m,n))/dx/(dx*12)) - ((BB(m,n))/12/dx) ;
        Am((m)*Nx,(n-1)*Nx+Nx) = -30*AA(m,n)/dx/(dx*12) + CC(m,n) ;
        Am((m)*Nx,(n-1)*Nx+Nx-1) = ((16*AA(m,n))/dx/(dx*12)) - (8*(BB(m,n))/12/dx);
        Am((m)*Nx,(n-1)*Nx+Nx-2) = (-(AA(m,n))/dx/(dx*12)) + ((BB(m,n))/12/dx);
    end

    A1 = As(:,:,2);
    B1 = Bs(:,:,2);
    C1 = Cs(:,:,2);
    AA = A1;
    BB = B1;
    CC = C1;

    for n=1:Nm
        Am((m-1)*Nx+2,(n-1)*Nx+1) = ((16*AA(m,n))/dx/(dx*12)) - (8*(BB(m,n))/12/dx);
        Am((m-1)*Nx+2,(n-1)*Nx+2) = (-30*(AA(m,n))/dx/(dx*12)) + (CC(m,n)) ;
        Am((m-1)*Nx+2,(n-1)*Nx+3) = ((16*AA(m,n))/dx/(dx*12)) + (8*(BB(m,n))/12/dx) ;
        Am((m-1)*Nx+2,(n-1)*Nx+4) = (-(AA(m,n))/dx/(dx*12)) - ((BB(m,n))/12/dx)  ;
        Am((m-1)*Nx+2,(n-1)*Nx+Nx-1) = (-(AA(m,n))/dx/(dx*12)) + ((BB(m,n))/12/dx);
    end
    A1 = As(:,:,Nx-1);
    B1 = Bs(:,:,Nx-1);
    C1 = Cs(:,:,Nx-1);
    AA = A1;
    BB = B1;
    CC = C1;

    for n=1:Nm
        Am((m)*Nx-1,(n-1)*Nx+2) = (-(AA(m,n))/dx/(dx*12)) - ((BB(m,n))/12/dx) ;
        Am((m)*Nx-1,(n-1)*Nx+Nx-3) = (-(AA(m,n))/dx/(dx*12)) + ((BB(m,n))/12/dx) ;
        Am((m)*Nx-1,(n-1)*Nx+Nx) = ((16*AA(m,n))/dx/(dx*12)) + (8*(BB(m,n))/12/dx) ;
        Am((m)*Nx-1,(n-1)*Nx+Nx-1) = (-30*(AA(m,n))/dx/(dx*12)) + (CC(m,n));
        Am((m)*Nx-1,(n-1)*Nx+Nx-2) = ((16*AA(m,n))/dx/(dx*12)) - (8*(BB(m,n))/12/dx);

    end



    for i1 = 3:Nx-2
        A = As(:,:,i1);B = Bs(:,:,i1); C = Cs(:,:,i1);

        AA = A;
        BB = B;
        CC = C;
        im=(m-1)*Nx+i1;
        for n=1:Nm

            i2 = i1-2;
            in = (n-1)*Nx+i2;
            Am(im,in) = (-(AA(m,n))/dx/(dx*12)) + ((BB(m,n))/12/dx);

            i2 = i1-1;
            in = (n-1)*Nx+i2;
            Am(im,in) = ((16*AA(m,n))/dx/(dx*12)) - (8*(BB(m,n))/12/dx) ;

            i2 = i1;
            in = (n-1)*Nx+i2;
            Am(im,in) = (-30*(AA(m,n))/dx/(dx*12)) + (CC(m,n));

            i2 = i1+1;
            in = (n-1)*Nx+i2;
            Am(im,in) = ((16*AA(m,n))/dx/(dx*12)) + (8*(BB(m,n))/12/dx)  ;

            i2 = i1+2;
            in = (n-1)*Nx+i2;
            Am(im,in) = (-(AA(m,n))/dx/(dx*12)) - ((BB(m,n))/12/dx)  ;


        end
    end
end

constr = zeros(Nx,Nm*Nx);
for i=1:Nm
    constr(1:end,(i-1)*(Nx)+1:i*Nx) = eye(Nx) ;
end

Bv((Nm-1)*Nx+1:end) = psi(1:Nx);
Am = sparse(Am);
constr = sparse(constr);

Amm = [Am;constr];
Xm = Amm\Bv;

ph = [Xm(1:Nx)];  % phi_{-2} modal coefficient
end



