function [r,z] = shape_edit(r0,z0,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [r,z] = shape_edit(r0,z0,s)
%
%  PURPOSE: Edit a plasma shape by specifying new shape parameters
%           Original shape parameters can be found with
%           s0 = shape_params(r0,z0)
%
%  INPUTS:  s, structure with new shape parameters, these can be edited:
%              rsurf, zsurf, aminor, bminor, elong, triu, tril,
%              squo, squi, sqli, sqlo
%           Omit fields or set value to nan for parameters to be ignored
%           One of aminor, bminor, elong must be ignored (default bminor)
%
%  OUTPUTS: r,z, modified boundary coordinates
%
%  See also SHAPE_CREATE, SHAPE_PARAMS
		
%
%  WRITTEN BY:  Anders Welander ON 2017-12-15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r,z will be modified
r = r0;
z = z0;

% Basis shape parameters and extras
b = shape_params(r,z);


% Edit basis shape
if isfield(s,'rsurf')
  r = r + s.rsurf - b.rsurf;
end
if isfield(s,'zsurf')
  z = z + s.zsurf - b.zsurf;
end
b = shape_params(r,z);
if isfield(s,'aminor')
  r = b.rsurf + (r-b.rsurf)*s.aminor/b.aminor;
end
if isfield(s,'bminor')
  z = b.zsurf + (z-b.zsurf)*s.bminor/b.bminor;
end
b = shape_params(r,z);
if isfield(s,'elong') 
  if ~isfield(s,'aminor')
    aminor = b.bminor/s.elong;
    r = b.rsurf + (r-b.rsurf)*aminor/b.aminor;
  else
    bminor = s.aminor*s.elong;
    z = b.zsurf + (z-b.zsurf)*bminor/b.bminor;
  end
end

% Update shape parameters in b and also get indices for triangularity modifications
[b,k1,k2,k3,k4] = shape_params(r,z);

% Limits on shifted points
rmin = b.R8(5)+b.aminor/1e3;
rmax = b.R8(1)-b.aminor/1e3;

if isfield(s,'triu')
  ru = min(rmax,max(rmin,b.rsurf - s.triu*b.aminor));
  f = (ru-b.R8(1))/(b.R8(3)-b.R8(1))-1;
  dr1 = r(k1)-b.R8(1);
  r(k1) = b.R8(1) + dr1.*(1+f*(z(k1)-b.Z8(1))/(b.Z8(3)-b.Z8(1)));
  f = (ru-b.R8(5))/(b.R8(3)-b.R8(5))-1;
  dr2 = r(k2)-b.R8(5);
  r(k2) = b.R8(5) + dr2.*(1+f*(z(k2)-b.Z8(5))/(b.Z8(3)-b.Z8(5)));
end
if isfield(s,'tril')
  rl = min(rmax,max(rmin,b.rsurf - s.tril*b.aminor));
  f = (rl-b.R8(5))/(b.R8(7)-b.R8(5))-1;
  dr3 = r(k3)-b.R8(5);
  r(k3) = b.R8(5) + dr3.*(1+f*(z(k3)-b.Z8(5))/(b.Z8(7)-b.Z8(5)));
  f = (rl-b.R8(1))/(b.R8(7)-b.R8(1))-1;
  dr4 = r(k4)-b.R8(1);
  r(k4) = b.R8(1) + dr4.*(1+f*(z(k4)-b.Z8(1))/(b.Z8(7)-b.Z8(1)));
end

% Update shape parameters in b and also get coordinates for squareness modifications
[b,k1,k2,k3,k4,x1,x2,x3,x4,y1,y2,y3,y4,i1,i2,i3,i4] = shape_params(r,z);

if isfield(s,'squo') & ~isnan(s.squo)
  
  % The two points that affect the calculation in shape_params:
  x = x1(i1-1:i1);
  y = y1(i1-1:i1);
  
  % Get sq which differs only slightly from s.squo
  % and fraction f of basis shape details to keep
  [sq,f] = tweaked_sq(x,y,b.squo,s.squo);
  
  % rho and angles for all points
  rho1_basis = sqrt(x1.^2+y1.^2);
  a1 = angle(complex(x1,y1));
  
  % For an ideal shape at the same squo, rho is
  rho1_ideal = ideal_shape(b.squo,a1);
  
  % The difference contains the details of the basis shape
  drho1 = rho1_basis - rho1_ideal;
  
  % The ideal target
  rho1_ideal_target = ideal_shape(sq,a1);
  
  % The new rho for all points
  rho1 = rho1_ideal_target + f*drho1;
  
  % Limit on rho1
  rho1(rho1 < 0) = 0;
  
  % New r(k1), z(k1)
  r1 = b.R8(3) + (b.R8(1)-b.R8(3))*rho1.*cos(a1);
  z1 = b.Z8(1) + (b.Z8(3)-b.Z8(1))*rho1.*sin(a1);
  
  % Limits on r,z
  r1(r1 > b.R8(1)-b.aminor/1e3) = b.R8(1)-b.aminor/1e3;
  z1(z1 > b.Z8(3)-b.bminor/1e3) = b.Z8(3)-b.bminor/1e3;
  
  % Restore first point which is shape point 1
  r1(1) = b.R8(1);
  z1(1) = b.Z8(1);
  
  % Update r, z
  r(k1) = r1;
  z(k1) = z1;
  
end

if isfield(s,'squi') & ~isnan(s.squi)
  
  % The two points that affect the calculation in shape_params:
  x = x2(i2-1:i2);
  y = y2(i2-1:i2);
  
  % Get sq which differs only slightly from s.squi
  % and fraction f of basis shape details to keep
  [sq,f] = tweaked_sq(x,y,b.squi,s.squi);
  
  % rho and angles for all points
  rho2_basis = sqrt(x2.^2+y2.^2);
  a2 = angle(complex(x2,y2));
  
  % For an ideal shape at the same squi, rho is
  rho2_ideal = ideal_shape(b.squi,a2);
  
  % The difference contains the details of the basis shape
  drho2 = rho2_basis - rho2_ideal;
  
  % The ideal target
  rho2_ideal_target = ideal_shape(sq,a2);
  
  % The new rho for all points
  rho2 = rho2_ideal_target + f*drho2;
  
  % Limit on rho2
  rho2(rho2 < 0) = 0;
  
  % New r(k2), z(k2)
  r2 = b.R8(3) + (b.R8(5)-b.R8(3))*rho2.*cos(a2);
  z2 = b.Z8(5) + (b.Z8(3)-b.Z8(5))*rho2.*sin(a2);
  
  % Limits on r,z
  r2(r2 < b.R8(5)+b.aminor/1e3) = b.R8(5)+b.aminor/1e3;
  z2(z2 > b.Z8(3)-b.bminor/1e3) = b.Z8(3)-b.bminor/1e3;
  
  % Restore first point which is shape point 3
  r2(1) = b.R8(3);
  z2(1) = b.Z8(3);
  
  % Update r, z
  r(k2) = r2;
  z(k2) = z2;
  
end

if isfield(s,'sqli') & ~isnan(s.sqli)
  
  % The two points that affect the calculation in shape_params:
  x = x3(i3-1:i3);
  y = y3(i3-1:i3);
  
  % Get sq which differs only slightly from s.sqli
  % and fraction f of basis shape details to keep
  [sq,f] = tweaked_sq(x,y,b.sqli,s.sqli);
  
  % rho and angles for all points
  rho3_basis = sqrt(x3.^2+y3.^2);
  a3 = angle(complex(x3,y3));
  
  % For an ideal shape at the same sqli, rho is
  rho3_ideal = ideal_shape(b.sqli,a3);
  
  % The difference contains the details of the basis shape
  drho3 = rho3_basis - rho3_ideal;
  
  % The ideal target
  rho3_ideal_target = ideal_shape(sq,a3);
  
  % The new rho for all points
  rho3 = rho3_ideal_target + f*drho3;
  
  % Limit on rho3
  rho3(rho3 < 0) = 0;
  
  % New r(k3), z(k3)
  r3 = b.R8(7) + (b.R8(5)-b.R8(7))*rho3.*cos(a3);
  z3 = b.Z8(5) + (b.Z8(7)-b.Z8(5))*rho3.*sin(a3);
  
  % Limits on r,z
  r3(r3 < b.R8(5)+b.aminor/1e3) = b.R8(5)+b.aminor/1e3;
  z3(z3 < b.Z8(7)+b.bminor/1e3) = b.Z8(7)+b.bminor/1e3;
  
  % Restore first point which is shape point 5
  r3(1) = b.R8(5);
  z3(1) = b.Z8(5);
  
  % Update r, z
  r(k3) = r3;
  z(k3) = z3;
  
end

if isfield(s,'sqlo') & ~isnan(s.sqlo)
  
  % The two points that affect the calculation in shape_params:
  x = x4(i4-1:i4);
  y = y4(i4-1:i4);
  
  % Get sq which differs only slightly from s.sqlo
  % and fraction f of basis shape details to keep
  [sq,f] = tweaked_sq(x,y,b.sqlo,s.sqlo);
  
  % rho and angles for all points
  rho4_basis = sqrt(x4.^2+y4.^2);
  a4 = angle(complex(x4,y4));
  
  % For an ideal shape at the same sqlo, rho is
  rho4_ideal = ideal_shape(b.sqlo,a4);
  
  % The difference contains the details of the basis shape
  drho4 = rho4_basis - rho4_ideal;
  
  % The ideal target
  rho4_ideal_target = ideal_shape(sq,a4);
  
  % The new rho for all points
  rho4 = rho4_ideal_target + f*drho4;
  
  % Limit on rho4
  rho4(rho4 < 0) = 0;
  
  % New r(k4), z(k4)
  r4 = b.R8(7) + (b.R8(1)-b.R8(7))*rho4.*cos(a4);
  z4 = b.Z8(1) + (b.Z8(7)-b.Z8(1))*rho4.*sin(a4);
  
  % Limits on r,z
  r4(r4 > b.R8(1)-b.aminor/1e3) = b.R8(1)-b.aminor/1e3;
  z4(z4 < b.Z8(7)+b.bminor/1e3) = b.Z8(7)+b.bminor/1e3;
  
  % Restore first point which is shape point 7
  r4(1) = b.R8(7);
  z4(1) = b.Z8(7);
  
  % Update r, z
  r(k4) = r4;
  z(k4) = z4;
  
end


% rho for ideal shape as function of squareness in square [0 1],[0 1]
% angles range from 0 to pi/2
function [rho, rhosquare] = ideal_shape(squareness,angles)

% Distance from circle to outer corner
CD = sqrt(2)-1;

% rho at 45 degrees, the point that determines squareness
rho45 = 1 + squareness*CD;

% Squareness = 1 is a square with rho45 = sqrt(2)
rhosquare = 1./cos(angles);
rhosquare(angles>pi/4) = 1./sin(angles(angles>pi/4));

% Squareness = -sqrt(2)/2 is a triangle with rho45 = sqrt(2)/2
rhotriangle = 1./(cos(angles)+sin(angles));

if rho45 < sqrt(2)/2
  % The shape has an indentation to achieve squareness below -sqrt(2)/2
  drhoindent = sqrt(2)/2 - rho45;
  f = drhoindent/(1-exp(-pi^2/128));
  rho = rhotriangle - f*exp(-(angles-pi/4).^2/8) + f*exp(-pi^2/128);
elseif rho45 < 1
  % Linear interpolation between the triangle and the circle
  f = (1-rho45)/(1-sqrt(2)/2);
  rho = f*rhotriangle + (1-f); % rhocircle = 1
elseif rho45 < 1.2
  % Formula for boundary bulging out from circle
  f = (rho45-1)/(1-exp(-pi^2/2));
  rho = 1 + f*exp(-8*(angles-pi/4).^2);
  rho(rho>rhosquare) = rhosquare(rho>rhosquare);
else
  % Linear interpolation for boundary morphing into a square
  f = 0.2/(1-exp(-pi^2/2));
  rho = 1 + f*exp(-8*(angles-pi/4).^2);
  rho(rho>rhosquare) = rhosquare(rho>rhosquare); % Solution for rho = 1.2
  f = (rho45-1.2)/(sqrt(2)-1.2); % Fraction of rhosquare in solution
  rho = f*rhosquare + (1-f)*rho; % A part rhosquare and a part rho(1.2)
end


% Finds small difference in sq between ideal_shape and shape_params
% x,y are two points on either side of 45 degrees
% sqbasis is squareness calculated with x,y
% sqtarget is desired new value of squareness
function [sq,f] = tweaked_sq(x,y,sqbasis,sqtarget)
  
  % rho and angles
  rho = sqrt(x.^2+y.^2);
  a = angle(complex(x,y));
  
  % For an ideal shape at the same squareness
  rho_ideal = ideal_shape(sqbasis,a);
  
  % The difference contains the details of the basis shape
  drho = rho - rho_ideal;
  
  % Fraction of drho to keep
  if abs(sqtarget) < 0.9
    if sqtarget > sqbasis
      f = 1 - (sqtarget-sqbasis)/(1-sqbasis); % Keep zero at squarenss = +1
    else
      f = 1 + (sqtarget-sqbasis)/(1+sqbasis); % Keep zero at squarenss = -1
    end
  else
    f = 1; % Keep all of drho
  end
  
  % Initial sq, will be slightly modified in loop below
  sq = sqtarget;
  
  % Iterate to remove a small error related to interpolation in shape_params
  for i = 1:33
  
    % The ideal target
    rho_ideal_target = ideal_shape(sq,a);

    % The values of x, y for the target shape
    rho_target = rho_ideal_target + f*drho;
    x = rho_target.*cos(a);
    y = rho_target.*sin(a);

    % The calculation of squareness that is done in shape_params
    g = (y(1)-x(1))/(x(2)-x(1)-y(2)+y(1));
    sqact = ((g*x(2)+(1-g)*x(1))*sqrt(2)-1)/(sqrt(2)-1);

    % Modify sq so that sqact more exactly equals sqtarget
    sq = sq + sqtarget-sqact;
    
    if abs(sqtarget-sqact) < 1e-14
      break
    end

  end 
  % Now sq is known

