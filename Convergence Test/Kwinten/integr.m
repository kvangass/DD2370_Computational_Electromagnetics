% --------------------------------------------------------------
% Compute potential on symmetry axis of square plate
% --------------------------------------------------------------

% Arguments:
%    z    = the height over the plate
%    a    = the side of the square
%    n    = the number of elements along each side of the plate
%    rule = a string 'midpoint' or 'simpson' that specifies
%           the integration rule
% Returns:
%    pot = the potential at the point (0,0,z)

function pot = integr(z, a, n, rule)

x  = linspace(0, a, n+1); 
y  = linspace(0, a, n+1);
h  = a/n;
zs = z^2;

if (strcmp(rule, 'midpoint'))
  profile on
  % Midpoint integration
  xs(1:n) = (x(1:n) + h/2).^2;
  ys(1:n) = (y(1:n) + h/2).^2;
  % Set the accumulator to zero
  intSquare = 0;
  for i=1:n
      for j=1:n
          intSquare = intSquare + 1./sqrt(xs(i) + ys(j) + zs); % integral of small square
      end
  end
  profile off
  

  
elseif (strcmp(rule, 'simpson'))
  profile on
  % Simpson's rule
  intSquare = 0;
  for i = 1:n
    x1 = x(i)^2; x2 = (x(i) + h/2)^2; x3 = (x(i) + h)^2;
    y1(1:n) = y(1:n).^2; 
    y2(1:n) = (y(1:n) + h/2).^2;
    y3(1:n) = (y(1:n) + h).^2;
    intSquare = intSquare + sum(  1./sqrt(x1+y1+zs) + 1./sqrt(x1+y3+zs) ...
                      + 1./sqrt(x3+y1+zs) + 1./sqrt(x3+y3+zs)...
                      + 4./sqrt(x2+y1+zs) + 4./sqrt(x2+y3+zs)...
                      + 4./sqrt(x1+y2+zs) + 4./sqrt(x3+y2+zs)...
                      + 16./sqrt(x2+y2+zs))/36;
  end
  profile off
else
  
  error(['Only midpoint integration and Simpson''s rule are ' ...
         'implemented'])
  
end
format long
pot = intSquare*h^2;

