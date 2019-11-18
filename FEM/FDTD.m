%Task 1
close all; clear all

N=200;
cap_fem = 92.7633;
cap_abs = capacitor_FD2(1,1,2,2,N,1E-3,1);

n_vals = 11:1:111;
%n_vals = [25,50,100,200];
cap_vals = zeros(length(n_vals),1);
i=0;
for n = n_vals
    i = i+1;
    cap_1 = capacitor_FD2(1,1,2,2,n,1E-3,1);
    cap_vals(i) = cap_1;
end

error = abs(cap_vals - cap_abs);
error_fem = abs(cap_fem - cap_abs);
error_fem = ones(length(n_vals),1)*error_fem;


figure(1)
semilogy(n_vals, error,'b')
hold on
semilogy(n_vals,error_fem,'r')
grid on
title('Error Comparison FDTD vs FEM')
ylabel('Error')
xlabel('N')
legend('FDTD','FEM')