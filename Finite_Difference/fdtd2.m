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
% Use Courant condition
dt=dx/(2.*c0);
% Constants
cc=c0*dt/dx;
% Initialize vectors
Ez=zeros(1,ke);
Hy=zeros(1,ke);
% Gaussian pulse
t0=20;
spread=8;
% Start loop
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

plot(Ez);axis([1 ke -2 2]);
xlabel('x');
ylabel('E_z');

pause()
end
