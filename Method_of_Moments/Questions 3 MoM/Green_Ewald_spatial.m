function Green = Green_Ewald_spatial( Dx,z,zp,k,kx,p,E,reg,eps )
%  Dx=x-x', z, z' are the spatial coordinates
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  reg            is 0 if the functions i NOT regularized, 1 if the function is regularized
%  eps            is the relative accuracy

ZERO= 1e-8;
err = 1;
Green = 0;


Rn = sqrt(Dx.^2 + (z-zp)^2);
aq=(0.5*k/E)^2;
x_eint=(Rn*E).^2;
if reg==0
    E_q= expint(x_eint);
    sum = E_q;
    factorBalanced=1;
    for q=1:41
        factorBalanced=factorBalanced*aq/q;
        E_qp1 = Next_Eq(E_q,q,x_eint);
        sum = sum + factorBalanced*E_qp1;
        E_q=E_qp1;
    end
    Green = sum;
else
    E_q= expint(x_eint);
    sum = E_q+log(x_eint);
    sum(Rn<ZERO) = -vpa(eulergamma);
    factorBalanced=1;
    for q=1:40
        factorBalanced=factorBalanced*aq/q;
        E_qp1 = Next_Eq(E_q,q,x_eint);
        E_qp1(Rn<ZERO) = 1/q;
        sum = sum + factorBalanced*E_qp1;
        E_q=E_qp1;
    end
    Green = sum;
end


n=1;
while max(err)>eps

    Rn = sqrt((Dx-n*p).^2 + (z-zp)^2);
    aq=(0.5*k/E)^2;
    x_eint=(Rn*E).^2;
    E_q= expint(x_eint);
    sum = E_q;
    factorBalanced=1;
    for q=1:40
        factorBalanced=factorBalanced*aq/q;
        E_qp1 = Next_Eq(E_q,q,x_eint);
        sum = sum + factorBalanced*E_qp1;
        E_q=E_qp1;
    end
    Green_temp = exp(-1i*n*kx*p)*sum;

    Rn = sqrt((Dx+n*p).^2 + (z-zp)^2);
    aq=(0.5*k/E)^2;
    x_eint=(Rn*E).^2;
    E_q= expint(x_eint);
    sum = E_q;
    factorBalanced=1;
    for q=1:40
        factorBalanced=factorBalanced*aq/q;
        E_qp1 = Next_Eq(E_q,q,x_eint);
        sum = sum + factorBalanced*E_qp1;
        E_q=E_qp1;
    end
    Green_temp = Green_temp + exp(1i*n*kx*p)*sum;

    n=n+1;

    Green = Green + Green_temp;
    err = abs((Green_temp)./Green);

end

Green = Green/(4*pi);


end

