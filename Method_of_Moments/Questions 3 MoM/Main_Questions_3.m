
clear all
close all

global k p kx eps reg E w N

mu0 = 4*pi*1e-7;
f  = 10e9;   % Frequency
k = 2*pi*f/3e8;
lambda = 2*pi/k;

method = 'PM';   % 'GA' for Galerkin, 'PM' for point matching

no_freq = 9;
Refl_coeff = zeros(1,no_freq);
for ind=1:no_freq

    ind    % Displays the step number

    p = ind*lambda/10;
    w = p/2;
    theta = 90*pi/180;
    kx = k*cos(theta);
    eps = 10e-8;
    reg = 0;
    
    E = sqrt(pi)/p;
    
    N = 20;    % Number of basis functions
    MoM_matrix = zeros(N);
    for n=1:N
        for m=1:N
            if (method=='PM')
                MoM_matrix(m,n) = MoM_element_PM( m,n );
            elseif (method=='GA')
                MoM_matrix(m,n) = MoM_element_Galerkin( m,n );
            else
                disp('Testing not available')
            end
        end
    end
    
    l=w/N;
    V_forcing = zeros(1,N);
    x_m = zeros(1,N);
    for m=1:N
        x_m(m) = l*(m-1/2)-w/2;
        if (method=='PM')
            V_forcing(m) = 0.0; %Write here the forcing term with Point Matching;
        elseif (method=='GA')
            V_forcing(m) = -l*exp(-1i*kx*x_m(m))*sinc(0.5*kx*l/pi)/(-1i*2*pi*f*mu0);
        else
            disp('Testing not available')
        end
    end
    
    I = mldivide(MoM_matrix,V_forcing.');    %equivalently   MoM_matrix\V_forcing
    
    % Computation of the reflection coefficient: Peterson, p.268 (7.39)
    reflection_term = zeros(1,N);
    for m=1:N
        reflection_term(m) = I(m)*exp(1i*kx*x_m(m));
    end
    Refl_coeff(ind) = -(2*pi*f*mu0*l)/(2*p*sqrt(k^2-kx^2))*sinc(0.5*kx*l/pi)*sum(reflection_term);



end

figure(1)
%Plot here the magnitude of the current on one strip
title('Magnitude of the current on one strip')

figure(2)
%Plot here the reflection coefficient
title('Reflection coefficient of the 0th harmonic')

