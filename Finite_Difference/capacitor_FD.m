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

% Make grids
h  = c/n;                % Grid size
na = round(0.5*a/h);         % Number of segments on 'a'
n;
m  = round(d/h);         % Number of segments on 'd'
h;
mb = round(0.5*b/h);         % Number of segments on 'b'

%x  = linspace(0,0.5*c,n+1);  % Grid points along x-axis
%y  = linspace(0,0.5*d,m+1);  % Grid points along y-axis

num_points = n*m;
A = zeros(num_points,num_points);
G = zeros(num_points,1);

%Generate mask for outer values then reshape
mask_outer = zeros(m,n);
mask_outer(1,:) = 1;
mask_outer(end,:) = 1;
mask_outer(:,1)=1;
mask_outer(:,end)=1;
flat_mask_outer = mask_outer(:);
%Generate mask for inner values then reshape
mask_inner = zeros(m,n);
m_center = round(m/2);
n_center = round(n/2);
mask_inner(m_center-mb:m_center+mb,n_center-na:n_center+na) = 1;
flat_mask_inner = mask_inner(:);
%Questions
%Q1: should we use the symmetry condition,
%Q2: should we use mirroring over the symmetry condition? j+1 -> j-1

%Itterate over all points in the grid
for row = 1:num_points
    %row is effectively the grid number in the f matrix
    if flat_mask_inner(row)%Inner Grid point
        A(row,row) = 1;
        G(row) = 1;
        
    elseif flat_mask_outer(row) %Outer grid point
        A(row,row) = 1;
        G(row) = 0;
        
    else %Poisson case
        G(row) = 0;
        %HERE , the diagonal element
        A(row,row) = -4/h^.2;

        %UP
        %i_up = i-1;
        %j_up = j;
        f_num_up = row -1;
        A(row,f_num_up) = 1/h.^2;

        %DOWN
        %i_down = i+1;
        %j_down = j;
        f_num_down = row + 1;
        A(row,f_num_down) = 1/h.^2;

        %Left
        %i_left = i;
        %j_left = j-1;
        f_num_left = row-m;
        A(row,f_num_left) = 1/h.^2;

        %Right
        %i_right = i;
        %j_right = j+1;
        f_num_right = row+m;
        A(row,f_num_right) = 1/h.^2;
    
    end
end

%Solve the linear system to find the values of phi
phi =A\G;
phi_grid = reshape(phi,m,n);

figure(1)
surf(phi_grid)

figure(2)
surf(mask_inner)

figure(3)
surf(mask_outer)
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
  q = q + (f(i,m)+f(i+1,m))*0.5; %integrate along upper boundary
end
for j = 1:m
  q = q + (f(n,j)+f(n,j+1))*0.5; %integrate along right boundary
end
cap = q*4;           % 4 quadrants
cap = cap*8.854187;  % epsilon0*1e12 gives answer in pF/m
