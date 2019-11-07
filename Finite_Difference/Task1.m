%Task 1
close all; clear all

n_vals = 5:5:100;
cap_vals = zeros(length(n_vals),1);
i=0;
for n = n_vals
    i = i+1;
    cap_1 = capacitor_FD2(0.5,0.5,2,2,n,1E-3,1);
    cap_vals(i) = cap_1;
end
%cap_2 = capacitor(0.5,0.5,2,2,30,1E-3,1)

figure(3)
plot(n_vals,cap_vals)
xlabel('n')
ylabel('cap')
title('Convergence')
