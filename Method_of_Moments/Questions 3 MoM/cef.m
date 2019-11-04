% Taken from http://dip.sun.ac.za/~weideman/research/cef.html

function w = cef(arg,N_pol)

%  Computes the function w(z) = exp(-z^2) erfc(-iz) using a rational 
%  series with N_pol terms.  It is assumed that Im(z) > 0 or Im(z) = 0.
%
%                             Andre Weideman, 1995

M = 2*N_pol;  M2 = 2*M;  k = [-M+1:1:M-1]';             % M2 = no. of sampling points.
L = sqrt(N_pol/sqrt(2));                                % Optimal choice of L.
theta = k*pi/M; t = L*tan(theta/2);                     % Define variables theta and t.
funct = exp(-t.^2).*(L^2+t.^2); funct = [0; funct];     % Function to be transformed.
a = real(fft(fftshift(funct)))/M2;                      % Coefficients of transform.
a = flipud(a(2:N_pol+1));                               % Reorder coefficients.
Z_arg = (L+1i*arg)./(L-1i*arg); pol = polyval(a,Z_arg); % Polynomial evaluation.
w = 2*pol./(L-1i*arg).^2+(1/sqrt(pi))./(L-1i*arg);      % Evaluate w(z).