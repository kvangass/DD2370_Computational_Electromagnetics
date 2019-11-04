function Green = Green_Ewald_spectral( Dx,z,zp,k,kx,p,E,eps )
%  Dx=x-x', z, z' are the spatial coordinates
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  E              Ewald parameter
%  eps            is the relative accuracy


err = 1;
kzn = sqrt(k^2-kx^2);
if imag(kzn)>0
    kzn=-kzn;
end
arg_erfc_p = 1i*kzn/(2*E)+abs(z-zp)*E;
erfc_p = exp(-arg_erfc_p^2)*cef(1i*arg_erfc_p,10);
arg_erfc_m = 1i*kzn/(2*E)-abs(z-zp)*E;
erfc_m = exp(-arg_erfc_m^2)*cef(1i*arg_erfc_m,10);
sum = erfc_p*exp(1i*kzn*abs(z-zp)) + erfc_m*exp(-1i*kzn*abs(z-zp));
Green = exp(-1i*kx*Dx)*sum/kzn;

n=1;
while err>eps
    
    kxn = kx + 2*pi*n/p;
    kzn = sqrt(k^2-kxn^2);
    if imag(kzn)>0
        kzn=-kzn;
    end
    arg_erfc_p = 1i*kzn/(2*E)+abs(z-zp)*E;
    erfc_p = exp(-arg_erfc_p^2)*cef(1i*arg_erfc_p,10);
    arg_erfc_m = 1i*kzn/(2*E)-abs(z-zp)*E;
    erfc_m = exp(-arg_erfc_m^2)*cef(1i*arg_erfc_m,10);
    sum = erfc_p*exp(1i*kzn*abs(z-zp)) + erfc_m*exp(-1i*kzn*abs(z-zp));
    Green_temp = exp(-1i*kxn*Dx)*sum/kzn;
    
    kxn = kx - 2*pi*n/p;
    kzn = sqrt(k^2-kxn^2);
    if imag(kzn)>0
        kzn=-kzn;
    end
    arg_erfc_p = 1i*kzn/(2*E)+abs(z-zp)*E;
    erfc_p = exp(-arg_erfc_p^2)*cef(1i*arg_erfc_p,10);
    arg_erfc_m = 1i*kzn/(2*E)-abs(z-zp)*E;
    erfc_m = exp(-arg_erfc_m^2)*cef(1i*arg_erfc_m,10);
    sum = erfc_p*exp(1i*kzn*abs(z-zp)) + erfc_m*exp(-1i*kzn*abs(z-zp));
    Green_temp = Green_temp + exp(-1i*kxn*Dx)*sum/kzn;
    
    n=n+1;
    
    Green = Green + Green_temp;
    err = abs((Green_temp)./Green);
    
end

Green = Green/(4*1i*p);


end