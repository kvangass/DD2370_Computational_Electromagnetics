%main

z=1;
a=5;
n=5;

pot1000 = integr(z, a, 1000, "simpson");

vals = [];
n_s = [20,40,60,80];
for n = n_s
    pot = integr(z, a, n, "simpson");
    vals = [vals pot];
end

error_Int = abs(vals - pot1000);

figure(1);
loglog(n_s, error_Int,'-');
hold on;
loglog(n_s, 1./n_s.^2 );
loglog(n_s, 1./n_s.^3 );
loglog(n_s, 1./n_s.^4 );
loglog(n_s, 1./n_s.^5 );

legend('Midpoint Error','2nd order', '3rd order','4th order','5th order')
xlabel("N");
ylabel("Error");
grid on;


h2 = 1./n_s.^2;
figure(2)
plot(h2,vals)
pfit = polyfit(h2, vals,2);
pfit
