function Green = Green_spectral( Dx,z,zp,k,kx,p,eps )
%  Dx=x-x', z, z' are the spatial coordinates
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  eps            is the relative accuracy

eps = 1e-8;
err = 1;
Green = 0;
kzn = sqrt(k^2-kx^2);
Green = exp(-1i*kx*Dx)*exp(-1i*kzn*abs(z-zp))/kzn;

% harmonics_to_plot = zeros(1,100);
% index = 1;

n=1;
while err>eps  %-1
    
    kxn = kx + 2*pi*n/p;
    kzn = sqrt(k^2-kxn^2);
    if imag(kzn)>0
        kzn=-kzn;
    end
    Green_temp = exp(-1i*kxn*Dx)*exp(-1i*kzn*abs(z-zp))/kzn;
%     harmonics_to_plot(index)= abs(Green_temp);
%     index = index+1;
    
    kxn = kx - 2*pi*n/p;
    kzn = sqrt(k^2-kxn^2);
    if imag(kzn)>0
        kzn=-kzn;
    end
    Green_temp = Green_temp + exp(-1i*kxn*Dx)*exp(-1i*kzn*abs(z-zp))/kzn;
    
    n=n+1;
    
    Green = Green + Green_temp;
    err = abs((Green_temp)./Green);
    
%     if(index>100)
%        plot(harmonics_to_plot) 
%        pause
%     end

    
end

Green = Green/(2*1i*p);




end