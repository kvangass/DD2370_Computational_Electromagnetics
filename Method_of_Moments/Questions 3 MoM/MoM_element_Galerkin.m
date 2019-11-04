function MoM_el = MoM_element_Galerkin( m,n )
%  n,m  are the MoM matrix element index
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  eps            is the relative accuracy

global E w N reg

l = w/N;
xm = l*(m-1/2)-w/2;
xn = l*(n-1/2)-w/2;
lim_inf = xm-l/2;
lim_sup = xm+l/2;
lim_inf_p = xn-l/2;
lim_sup_p = xn+l/2;

if n==m
    reg=1;
    MoM_el = integral2(@Green_function_Galerkin,lim_inf,lim_sup,lim_inf_p,lim_sup_p);
    MoM_el = MoM_el - (l^2/(2*pi))*(log(l*E)-3/2);
else
    reg=0;
    MoM_el = integral2(@Green_function_Galerkin,lim_inf,lim_sup,lim_inf_p,lim_sup_p);
end


end

