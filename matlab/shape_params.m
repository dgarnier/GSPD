function [s,k1,k2,k3,k4,x1,x2,x3,x4,y1,y2,y3,y4,i1,i2,i3,i4] = shape_params(r,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   s = shape_params(r,z)
%
%  PURPOSE: Return shape parameters for a boundary
%
%  INPUTS:  r, z, boundary coordinates
%
%  OUTPUTS: s, structure with shape parameters
%
%  See also SHAPE_CREATE, SHAPE_EDIT
		
%
%  WRITTEN BY:  Anders Welander ON 2017-12-15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of boundary points
n = length(r);

% Area & boundary direction
A = (sum((r(1:n-1)+r(2:n)).*(z(2:n)-z(1:n-1))));
dir = sign(A);
A = abs(A);
if dir == 0
  d.R8 = 'R for boundary points that define shape parameters';
  s.R8 = nan;
  d.Z8 = 'Z for boundary points that define shape parameters';
  s.Z8 = nan;
  d.rsurf = 'R for geometric center';
  s.rsurf = nan;
  d.zsurf = 'Z for geometric center';
  s.zsurf = nan;
  d.aminor = 'Half the radial width of the plasma';
  s.aminor = nan;
  d.bminor = 'Half the vertical height of the plasma';
  s.bminor = nan;
  d.elong = 'Elongation = bminor/aminor';
  s.elong = nan;
  d.triu = 'Upper triangularity';
  s.triu = nan;
  d.tril = 'Lower triangularity';
  s.tril = nan;
  d.squo = 'Upper outer squareness';
  s.squo = nan;
  d.squi = 'Upper inner squareness';
  s.squi = nan;
  d.sqli = 'Lower inner squareness';
  s.sqli = nan;
  d.sqlo = 'Lower outer squareness';
  s.sqlo = nan;
  %  error('r,z is not a boundary')
  return
end

% Find extremes
[~,io] = max(r);
[~,iu] = max(z);
[~,ii] = min(r);
[~,il] = min(z);

% Indices in each quadrant
if dir == 1 % indices go in positive theta direction
  if iu > io
    k1 = [io:iu-1];
  else
    k1 = [io:n 1:iu-1];
  end
  if ii > iu
    k2 = [iu:ii-1];
  else
    k2 = [iu:n 1:ii-1];
  end
  if il > ii
    k3 = [ii:il-1];
  else
    k3 = [ii:n 1:il-1];
  end
  if io > il
    k4 = [il:io-1];
  else
    k4 = [il:n 1:io-1];
  end
else % indices go in negative theta direction
  if iu < io
    k1 = [io:-1:iu+1];
  else
    k1 = [io:-1:1 n:-1:iu+1];
  end
  if ii < iu
    k2 = [iu:-1:ii+1];
  else
    k2 = [iu:-1:1 n:-1:ii+1];
  end
  if il < ii
    k3 = [ii:-1:il+1];
  else
    k3 = [ii:-1:1 n:-1:il+1];
  end
  if io < il
    k4 = [il:-1:io+1];
  else
    k4 = [il:-1:1 n:-1:io+1];
  end
end

% Fill in half the shape points
R8 = [r(io) nan r(iu) nan r(ii) nan r(il) nan];
Z8 = [z(io) nan z(iu) nan z(ii) nan z(il) nan];

% Position
rsurf = (R8(1)+R8(5))/2;
zsurf = (Z8(3)+Z8(7))/2;

% Size
aminor = (R8(1)-R8(5))/2;
bminor = (Z8(3)-Z8(7))/2;

% Elongation
elong = bminor/aminor;

% Triangularity
triu = (rsurf-R8(3))/aminor;
tril = (rsurf-R8(7))/aminor;

% Squareness
x1 = (r(k1)-R8(3))/(R8(1)-R8(3));
y1 = (z(k1)-Z8(1))/(Z8(3)-Z8(1));
for i1 = 2:numel(x1) % Starting at shape point 1, going to 3
  if y1(i1) > x1(i1)
    break % passed the point at 45 degrees
  end
end
% f1*x1(i1)+(1-f1)*x1(i1-1) = f1*y1(i1)+(1-f1)*y1(i1-1)
f1 = (y1(i1-1)-x1(i1-1))/(x1(i1)-x1(i1-1)-y1(i1)+y1(i1-1));
b1 = f1*x1(i1)+(1-f1)*x1(i1-1);
R8(2) = b1*(R8(1)-R8(3)) + R8(3);
Z8(2) = b1*(Z8(3)-Z8(1)) + Z8(1);
squo = (b1*sqrt(2)-1)/(sqrt(2)-1);

x2 = (r(k2)-R8(3))/(R8(5)-R8(3));
y2 = (z(k2)-Z8(5))/(Z8(3)-Z8(5));
for i2 = 2:numel(x2) % Starting at shape point 3, going to 5
  if x2(i2) > y2(i2)
    break % passed the point at 135 degrees
  end
end
f2 = (y2(i2-1)-x2(i2-1))/(x2(i2)-x2(i2-1)-y2(i2)+y2(i2-1));
b2 = f2*x2(i2)+(1-f2)*x2(i2-1);
R8(4) = b2*(R8(5)-R8(3)) + R8(3);
Z8(4) = b2*(Z8(3)-Z8(5)) + Z8(5);
squi = (b2*sqrt(2)-1)/(sqrt(2)-1);

x3 = (r(k3)-R8(7))/(R8(5)-R8(7));
y3 = (z(k3)-Z8(5))/(Z8(7)-Z8(5));
for i3 = 2:numel(x3) % Starting at shape point 5, going to 7
  if y3(i3) > x3(i3)
    break % passed the point at 225 degrees
  end
end
f3 = (y3(i3-1)-x3(i3-1))/(x3(i3)-x3(i3-1)-y3(i3)+y3(i3-1));
b3 = f3*x3(i3)+(1-f3)*x3(i3-1);
R8(6) = b3*(R8(5)-R8(7)) + R8(7);
Z8(6) = b3*(Z8(7)-Z8(5)) + Z8(5);
sqli = (b3*sqrt(2)-1)/(sqrt(2)-1);

x4 = (r(k4)-R8(7))/(R8(1)-R8(7));
y4 = (z(k4)-Z8(1))/(Z8(7)-Z8(1));
for i4 = 2:numel(x4) % Starting at shape point 7, going to 1
  if x4(i4) > y4(i4)
    break % passed the point at 315 degrees
  end
end
f4 = (y4(i4-1)-x4(i4-1))/(x4(i4)-x4(i4-1)-y4(i4)+y4(i4-1));
b4 = f4*x4(i4)+(1-f4)*x4(i4-1);
R8(8) = b4*(R8(1)-R8(7)) + R8(7);
Z8(8) = b4*(Z8(7)-Z8(1)) + Z8(1);
sqlo = (b4*sqrt(2)-1)/(sqrt(2)-1);


% Archive

d.R8 = 'R for boundary points that define shape parameters';
s.R8 = R8;

d.Z8 = 'Z for boundary points that define shape parameters';
s.Z8 = Z8;

d.rsurf = 'R for geometric center';
s.rsurf = rsurf;

d.zsurf = 'Z for geometric center';
s.zsurf = zsurf;

d.aminor = 'Half the radial width of the plasma';
s.aminor = aminor;

d.bminor = 'Half the vertical height of the plasma';
s.bminor = bminor;

d.elong = 'Elongation = bminor/aminor';
s.elong = elong;

d.triu = 'Upper triangularity';
s.triu = triu;

d.tril = 'Lower triangularity';
s.tril = tril;

d.squo = 'Upper outer squareness';
s.squo = squo;

d.squi = 'Upper inner squareness';
s.squi = squi;

d.sqli = 'Lower inner squareness';
s.sqli = sqli;

d.sqlo = 'Lower outer squareness';
s.sqlo = sqlo;

s.info = d;

