function [ As, Bs, Cs ,k ] = coef_ABC_flat( etas , h , Nm , mi0 ,dx ,h0  )
Nx = length(etas);
etasx = gradientp(etas,dx);
etasx2 = gradientp(etasx,dx);
hx = gradientp(h,dx);
hx2 = gradientp(hx,dx);

H = etas + h;
Hx = etasx + hx;
Hxx = etasx2 + hx2;




% Calculation of local wavenumbers 

k = zeros(Nx,Nm-2);
mi = mi0*H;
for i = 1:Nx

    % Propagating k0(x)

    % Initial guess using the Padé scheme from Hunt 1979
    kh0 = mi(i);
    d1 = 0.6666666666;
    d2 = 0.3555555555;
    d3 = 0.1608465608;
    d4 = 0.0632098765;
    d5 = 0.0217540484;
    d6 = 0.0065407983;
    den = 1 + d1*kh0 + d2*kh0.^2 + d3*kh0.^3 + d4*kh0.^4 + d5*kh0.^5 + d6*kh0.^6;
    k00 = sqrt((kh0.^2)+(kh0./den));
    x0 = k00*H(i);
    % Newton-Raphson
    for kk=1:300
        x = x0 - (x0*tanh(x0)-mi(i))/(tanh(x0) + x0*sech(x0)^2);
        if(abs((x - x0)/x) < 10^(-15))              %If the result is within the desired tolerance
        break;                                        %Done, so leave the loop
        end
    x0 = x;
    end
    kn(1) = x0/H(i);
    k(i,1) = kn(1);
    
    % Evanescent wavenumbers
    for n = 2:Nm-2

        % Initial guess (2nd order semi-explicit)
        kn0 = (n-1)*pi - ...
                ( (mi(i)^2 + ((n-1)*pi)^2 - mi(i))/(mi(i)^2 + ((n-1)*pi)^2 - 2*mi(i)))*atan(mi(i)/(n-1)/pi);    
        
        % Newton-Raphson
        x0 = kn0;
        for kk=1:500
            if (abs(tan(x0)+x0*sec(x0)^2)<1*10^(-15))
                x=x0;
                break;
            end
            x = x0 - (x0*tan(x0)+mi(i))/(tan(x0)+x0*sec(x0)^2);
            if(abs(x - x0)/abs(x) < 1*10^(-15))              
                break;                                       
            end
            x0 = x;
        end
        kn(n) = x0/H(i);
    end
    k(i,:) = kn;
end

% Analytical calculation of the spatial derivatives of local wavenumbers
k0 = abs(k(:,1));
DHk = - k0.*(k0.^2 - mi0^2)./(mi0 + H.*(k0.^2 - mi0^2));
kx(:,1) = DHk.*Hx;

DHHk = (k0.*(k0.^2 - mi0^2).^2)./(mi0 + H.*(k0.^2 - mi0^2)).^2  +...
    DHk.*(-3*(k0.^2)*mi0 + mi0.^3 - H.*(k0.^2 - mi0.^2).^2)./(mi0 + H.*(k0.^2 - mi0^2)).^2;

kxx(:,1) = DHHk.*Hx.^2 + DHk.*Hxx;

for n = 4:Nm
    DHk = -(k(:,n-2).*(k(:,n-2).^2 + mi0^2)./ (-mi0 + H.*(k(:,n-2).^2 + mi0^2)) );
    kx(:,n-2) = DHk.*Hx;
    kev = k(:,n-2);
    kxx(:,n-2) = -((2*(kev.^2).*Hx.*kx(:,n-2))./(-mi0 + H.*(mi0^2 + kev.^2))) -...
        ((mi0^2 + kev.^2).*Hx.*kx(:,n-2))./(-mi0 + H.*(mi0^2 + kev.^2)) +...
        (kev.*(mi0^2 + kev.^2).*Hx.*((mi0^2 + kev.^2).*Hx + 2*H.*kev.*kx(:,n-2)))./(-mi0 + H.*(mi0^2 + kev.^2)).^2 -...
        kev.*(mi0^2 + kev.^2).*Hxx./(-mi0 + H.*(mi0^2 + kev.^2));
end


As =  Acoefs(etas , h  ,mi0, h0 , k , Nm);
Bs =  Bcoefs(etas , h , etasx , hx ,mi0, h0 , k , kx , As, Nm);
Cs =  Ccoefs(etas , h , etasx , hx ,etasx2, hx2 ,mi0, h0 , k , kx ,kxx, As, Nm);
As(2,:,:)=[];
As(:,2,:)=[];
Bs(2,:,:)=[];
Bs(:,2,:)=[];
Cs(2,:,:)=[];
Cs(:,2,:)=[];

  
end


