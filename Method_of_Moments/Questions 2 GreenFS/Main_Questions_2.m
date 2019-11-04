clear all
close all

global k p kx eps reg E w N


% General parameters
mu0 = 4*pi*1e-7;    % Free-space permeability
f  = 10e9;   % Frequency
k = 2*pi*f/3e8;   % Free-space wavenumber
lambda = 2*pi/k;  % Free-space wavelength
reg = 0;          % Do not change

p = lambda/3; %lambda/10;    % Spatial period
theta = 30*pi/180;           % Incidence angle in radians (see Picture in slides)
kx = k*cos(theta);           % Wavenumber along the x direction (see Picture in slides)
eps = 10e-8;                 % Accuracy
    


%% Question 2.1
z = lambda/5; %0.0;
zp = 0.0;
Dx=-0.5*p:p/100:0.5*p;   % Dx=x'xp
Green_samples_spectral = zeros(1,length(Dx));
for index=1:length(Dx)
   % Computation of the spectral series
   Green_samples_spectral(index) = Green_spectral( Dx(index),z,zp,k,kx,p,eps );
end
figure(1)
plot(Dx,real(Green_samples_spectral))
title('Real part of the Green function (spectral series) in -0.5*p,0.5*p')
figure(2)
plot(Dx,imag(Green_samples_spectral))
title('Imaginary part of the Green function (spectral series) in -0.5*p,0.5*p')


%% Questions 2.2 and 2.3
%z = 0.0;          % Choose one value of z
z = lambda/5;    % Choose one value of z
zp = 0.0;
Dx= 0.5*p;   % Dx=x'xp
Green_spectral_harmonic_decay( Dx,z,zp,k,kx,p );


%% Question 2.4

% Setting the parameter
z = lambda/5; %0.0;
zp = 0.0;
Dx=-0.5*p:p/100:0.5*p;   % Dx=x'xp

% Computing the spectral series
Green_samples_spectral = zeros(1,length(Dx));
for index=1:length(Dx)
   % Computation of the spectral series
   Green_samples_spectral(index) = Green_spectral( Dx(index),z,zp,k,kx,p,eps );
end

% Computing the Ewald series
E = sqrt(pi)/p;
Green_samples_Ewald = zeros(1,length(Dx));
for index=1:length(Dx)
    % Computation of the Ewald series
   Green_samples_Ewald(index) = Green_Ewald_spatial( Dx(index),z,zp,k,kx,p,E,reg,eps ) + Green_Ewald_spectral( Dx(index),z,zp,k,kx,p,E,eps );
end

% Comparing the results
figure(4)
plot(Dx,real(Green_samples_Ewald),'r')
hold on
plot(Dx,real(Green_samples_spectral),'b')
title('Comparison spectral vs Ewald (Real parts in -0.5*p,0.5*p)')
figure(5)
plot(Dx,imag(Green_samples_Ewald),'r')
hold on
plot(Dx,imag(Green_samples_spectral),'b')
title('Comparison spectral vs Ewald (Imaginary parts in -0.5*p,0.5*p)')



%% Question 2.5

% Setting the parameter
z = 0.0;
zp = 0.0;
Dx=-0.5*p:p/100:0.5*p;   % Dx=x'xp
% Computing the Green's function
Green_samples_test = zeros(1,length(Dx));
for index=1:length(Dx)
    % Computation of the Green's function: which series do you use?!?
    Green_samples_test(index) = 0.0; %Replace the 0 with the series you want to use
end

% Comparing the results
figure(6)
plot(Dx,real(Green_samples_test))
title('Real part of the Green function for z=zp=0 and x in -0.5*p,0.5*p')
figure(7)
plot(Dx,imag(Green_samples_test))
title('Imaginary part of the Green function for z=zp=0 and x in -0.5*p,0.5*p')

