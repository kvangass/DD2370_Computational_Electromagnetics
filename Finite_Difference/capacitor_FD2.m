% --------------------------------------------------------------
% Compute capacitance per unit length of 
% a coaxial pair of rectangles
% --------------------------------------------------------------
function cap = capacitor_FD2(a, b, c, d, n, tol, rel)

% Arguments:
%    a   =  width of inner conductor
%    b   = height of inner conductor
%    c   =  width of outer conductor
%    d   = height of outer conductor
%    n   = number of points in the x-direction (horizontal)
%    tol = relative tolerance for capacitance
%    rel = relaxation parameter 
%          (optimum is 2-c/n, where c is about pi)  
% Returns:
%    cap = capacitance per unit length [pF/m]
% %%%%%%%%%%%%%%%%%%%%%
%Make Grids
N = n+1;
%Grid spacing
h = 0.5*c/n;
%number of segments on 'd'
m = round(0.5*d/h);

% N = n+1;
% n2 = a*n/c;
% n2=round(n2);
% h = a/n2;
% %number of segments on 'd'
% m = round(0.5*d/h);

%%%%%%%%%%%%%%%%%%%%%%%
M = m+1;
%Grid points along x-axis
x = linspace(0,0.5*c,N);
%Grid points along y-axis
y = linspace(0,0.5*d,M);
PhiCentral = 1; %Phi on the inner conductor

InConductor  = @(idk) ((d-b)/2 <= mod(idk-1,M)*h) && (a/2 >= h*floor(idk/(M+1)));

%make g in Ax = g
g = zeros((N)*(M),1);

for idk = 1:length(g)
    if (InConductor(idk))
        g(idk) = PhiCentral;
    end
end

Dx2 = sparse(length(g),length(g));
for idk = 1:length(g)
    [xi,yi] = VecToMatrix(idk,M);
    if(yi==1 || xi ==N) %handles top and right boundary condition f
        Dx2(idk,idk) = 1;
    elseif (xi == 1 ) %Handles left symmetry boundary
        Dx2(idk,idk) = -2;
        Dx2(idk,idk+M) = 2;
    elseif (InConductor(idk)) %handles inner conductor
        Dx2(idk,idk) = h^2/2; %(Dx2 + Dy2)/h^2 = 1 for this case
    else
        Dx2(idk,idk) = -2;
        Dx2(idk,idk+M) = 1;
        Dx2(idk,idk-M) = 1;
    end
end

Dy2 = sparse(length(g),length(g));
for idk = 1:length(g)
    [xi,yi] = VecToMatrix(idk,M);
    if (yi==1 || xi==N) %Handels top and right boundary
        Dy2(idk,idk) = 1;
    elseif (yi == M) %Handles the bottom symmetry
        Dy2(idk,idk) = -2;
        Dy2(idk,idk-1)= 2;
    elseif (InConductor(idk))
        Dy2(idk,idk) = h^2/2;
    else
        Dy2(idk,idk) = -2;
        Dy2(idk,idk-1) = 1;
        Dy2(idk,idk+1) = 1;
    end
end

Dxy2 = (Dx2+Dy2)/(h^2);
phi = Dxy2\g;

phi_2d = reshape(phi,M,N);
%figure(2)
%surf(phi_2d)

cap = gauss(n,m,h,phi_2d);
%End Capacitor Function
end

function [xi,yi] = VecToMatrix(idk,M)
    yi = mod(idk - 1,M) + 1;
    xi = ceil(idk/M);
end

% --------------------------------------------------------------
% Compute capacitance from the potential
% --------------------------------------------------------------
function cap = gauss(n,m,h,f)

% Arguments:
%    n    = number of points in the x-direction (horizontal)
%    m    = number of points in the y-direction (vertical)
%    h    = cell size
%    f    = 2D-array with solution
% Returns:
%    cap = capacitance per unit length [pF/m]
q = 0;
for i = 1:n
  q = q + (f(2,i)+f(2,i+1))*0.5; %integrate along upper boundary
end
for j = 2:m %Extra point?
  q = q + (f(j,n)+f(j+1,n))*0.5; %integrate along right boundary
end
cap = q*4;           % 4 quadrants
cap = cap*8.854187;  % epsilon0*1e12 gives answer in pF/m
end
