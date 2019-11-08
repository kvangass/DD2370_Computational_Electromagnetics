%Task 1
close all; clear all

cap_abs = capacitor_FD2(0.5,0.5,2,2,200,1E-3,1);

n_vals = 20:4:80;
%n_vals = [25,50,100,200];
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


error = abs(cap_vals - cap_abs);
figure(4)
loglog(n_vals,1./n_vals.^1,'r');
hold on;
loglog(n_vals,1./n_vals.^2,'g');
loglog(n_vals,1./n_vals.^3,'b');
loglog(n_vals,error,'k')
hold off;
grid on;
xlabel('N')
ylabel('Error')
