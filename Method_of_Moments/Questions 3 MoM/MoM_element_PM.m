function MoM_el = MoM_element_PM( m,n )
%  n,m  are the MoM matrix element index
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  eps            is the relative accuracy

global E w N reg xm

l = w/N;
xm = l*(m-1/2)-w/2;
xn = l*(n-1/2)-w/2;
lim_inf_p = xn-l/2;
lim_sup_p = xn+l/2;

reg=0;   % No extraction is necessary

%MoM_el = quad(@Green_function_PM,lim_inf_p,lim_sup_p);     % Numerical integration for continuous functions 

MoM_el = integral(@Green_function_PM,lim_inf_p,lim_sup_p);  % Numerical integration for possibly singular functions 

end

