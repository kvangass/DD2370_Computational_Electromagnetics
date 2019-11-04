function Green = Green_function_Galerkin( x,xp )
%  x, x', z, z'   are the spatial coordinates
%  k              is the medium wavenumber
%  kx             is the phase shift
%  p              is the spatial period
%  reg            is 0 if the functions i NOT regularized, 1 if the function is regularized
%  eps            is the relative accuracy

global k p kx eps reg E

z  = 0;
zp = 0;

Green = Green_Ewald_spatial( x-xp,z,zp,k,kx,p,E,reg,eps ) + Green_Ewald_spectral( x-xp,z,zp,k,kx,p,E,eps );



end

