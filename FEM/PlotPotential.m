% --------------------------------------------------------------
% Plot a 3D image of the potential distribution in the coaxial C.S.
% --------------------------------------------------------------
function PlotPotential(no2xy, el2no, sol)
% Arguments:
% no2xy = x- and y-coordinates of the nodes
% el2no = node indices for all triangles
% el2ed = edge indices for all elements
% sol   = Potential value
% Returns:
%-

% Sort the nodes of each element
el2no = sort(el2no);
ax = newplot;
no2xyz = [no2xy;sol'];
figure(1)
 patch('faces',el2no','vertices',no2xyz','facevertexcdata',sol, ...
 'facecolor',get(ax,'defaultsurfacefacecolor'), ...
 'edgecolor',get(ax,'defaultsurfaceedgecolor'));
hold on
  view(3)
xlabel('x [m]')
ylabel('y [m]')
zlabel('\phi [V]')
set(gca,'FontSize',15)
grid ON

% Plot the field itself as arrows
    
% Plot the 2D mesh
xy1 = no2xy(:,el2no(1,:));
xy2 = no2xy(:,el2no(2,:));
xy3 = no2xy(:,el2no(3,:));
xy = [xy1; xy2; xy3; xy1; NaN*xy1];
figure(2)
plot(xy(1:2:end),xy(2:2:end),'k')

no2xyz = [no2xy;sol'];
xyz1 = no2xyz(:,el2no(1,:));
xyz2 = no2xyz(:,el2no(2,:));
xyz3 = no2xyz(:,el2no(3,:));
xyz = [xyz1; xyz2; xyz3; xyz1; NaN*xyz1];

%Plot the 
figure(3)
plot3(xyz(1:3:end),xyz(2:3:end), xyz(3:3:end),'k')


 end