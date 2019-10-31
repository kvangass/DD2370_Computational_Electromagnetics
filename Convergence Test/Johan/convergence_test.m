%% Convergence
start=tic();
N = 73;
pot = zeros(1,N);
h = zeros(1,N);
for i=1:N
    pot(i) = integr(1, 1, 17 + i, 'simpson');
    h(i) = 1/(17 + i);
end
plot(h.^4,pot,'o')
grid ON
xlabel('h^{4}')
ylabel('\Phi')
time=toc(start)
%% error test

pot1000 = integr(1, 1, 1000, 'simpson');

errorInt = pot-pot1000;

loglog(1./h, errorInt)
hold on
loglog([20 40 60 80],1./[20 40 60 80].^2)
loglog([20 40 60 80],1./[20 40 60 80].^3)
loglog([20 40 60 80],1./[20 40 60 80].^4)
loglog([20 40 60 80],1./[20 40 60 80].^5)
xlabel('N')
ylabel('Error')
grid ON
legend('Simpson Error', '2nd order', '3rd order', '4th order','5th order')
%% Extrapolation
pfit = polyfit(h.^4,pot,3)
