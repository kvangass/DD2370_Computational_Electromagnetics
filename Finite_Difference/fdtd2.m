% This is a 1D FDTD simulation with pulse
% It displays a "movie" of the signal
% Size of the FDTD space
clear;
ke=50;
% Position of the source
ks=ke/2;
% Number of time steps
nsteps=50;
% Cell size and time stepping
c0=3.e8;
dx=0.01;
% Use Courant cond  ition
dt=dx/(2*c0);
% Constants
cc=c0*dt/dx;



%% E- and H-fields
% Initialize vectors
Ez=zeros(1,ke);
Hy=zeros(1,ke);
% Gaussian pulse
t0=20;
spread=8;
% Start loop
ExHist_EH = [];
for t=1:nsteps
   % E field loop
   for k=2:ke-1
    Ez(k)=Ez(k)+cc*(Hy(k-1)-Hy(k));
   end
   % Source
   Ez(ks)=exp(-.5*((t-t0)/spread)^2);
   % H field loop
   for k=1:ke-1
      Hy(k)=Hy(k)+cc*(Ez(k)-Ez(k+1));
   end
   ExHist_EH = [ExHist_EH; Ez];
plot(Ez);axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');

pause(0.1)
end

%% E-field using curl curl
% Initialize vectors
Ez=zeros(1,ke);
Ez1=zeros(1,ke);
Hy=zeros(1,ke);
Ez0 = Ez;
% Gaussian pulse
t0=20;
spread=8;
% Start loop
ExHist_E = [];
for t=1:nsteps
   % E field loop
   Ez_temp = Ez;
   for k=2:ke-1
    Ez1(k) = 2*Ez(k) - Ez0(k)  + cc^2*(Ez(k+1) - 2*Ez(k) + Ez(k-1)); 
   end
   
   % Source
   Ez0 = Ez_temp;
   Ez1(ks)=exp(-.5*((t-t0)/spread)^2);
   Ez = Ez1;
   ExHist_E = [ExHist_E; Ez];
plot(Ez);axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');

pause(0.1)
end


%% plot pulse propagation and dispersion

Fs_x = 2*pi/dx;
dF_x = Fs_x/ke;
kx = 0:dF_x:Fs_x/2-dF_x;

Fs_x = 1/dt;
dF_x = Fs_x/nsteps;
f = 0:dF_x:Fs_x/2-dF_x;
omega =  2*pi*f;


figure(1)
spectrum = fft2(ExHist_EH);
plot_spectrum = log(abs(spectrum));
pcolor(kx*dx, omega*dx/c0, plot_spectrum(1:ke/2,1:nsteps/2))
shading interp
title('Using E- and H-field')
xlabel('kx\Deltax')
ylabel('\omega\Deltax/c')
set(gca,'FontSize',15)

figure(2)
spectrum = fft2(ExHist_E);
plot_spectrum = log(abs(spectrum));
pcolor(kx*dx, omega*dx/c0, plot_spectrum(1:ke/2,1:nsteps/2))
shading interp
title('Using only E-field')
xlabel('kx\Deltax')
ylabel('\omega\Deltax/c')
set(gca,'FontSize',15)


figure(3)
plot(ExHist_E(12,:),'LineWidth',1.2);
hold on
plot(ExHist_E(24,:),'LineWidth',1.2);
plot(ExHist_E(36,:),'LineWidth',1.2);
plot(ExHist_E(48,:),'LineWidth',1.2);
axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');
title('Using only E-field')
legend('12 tstep','24 tstep', '36 tstep', '48 tstep')
set(gca,'FontSize',15)
hold off

figure(4)
plot(ExHist_EH(12,:),'LineWidth',1.2);
hold on
plot(ExHist_EH(24,:),'LineWidth',1.2);
plot(ExHist_EH(36,:),'LineWidth',1.2);
plot(ExHist_EH(48,:),'LineWidth',1.2);
axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');
title('Using E- and H-field')
legend('12 tstep','24 tstep', '36 tstep', '48 tstep')
set(gca,'FontSize',15)
hold off