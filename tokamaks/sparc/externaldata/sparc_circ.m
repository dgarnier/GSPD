% Defines various circuit connections, groupings, index labels, transition 
% matrices, etc. 
%
% This file is sparc-specific and will need to be modified if circuit
% connections are modified (for example, if vacuum vessel elements are
% grouped differently). 
% 
% Many of these circuit quantities are defined elsewhere in toksys. This is 
% an attempt to consolidate much of that info into one object. 
% Form the grouping circuits for the coils and vacuum vessel elements. 
% 
% cccirc is coil circuit connections. For example: [1 -1 2 2 2 3] says that
% coils 1 and 2 are connected in antiseries, coils 3-5 are connected in
% series, and coil 6 is independent. 
%
% vvgroup and vvcirc are similar, but subtly different. vvgroup = [1 1 2 3]
% says that elements 1 and 2 are lumped together (in parallel) into
% group 1, element 3 = group 2, and element 4 = group 3. vvcirc = [1 2 2]
% now applies on the groups, and says that group 1 is independent and
% groups 2 and 3 are connected in series. However, there are almost never 
% series connections in the vessel groupings, implying that for the majority
% of cases, vvcirc = 1:max(vvgroup); 
%
% Pcc, Pvv, and Pxx are transition maps from connected to unconnected
% circuits. Let ic, iv be unconnected coil and vessel current vectors. 
% Let icx, ivx be connected coil circuit and vessel circuit vectors. Let 
% ix := [ic; iv; ip] with ip=plasma current, ixx = [icx; ivx; ip]. 
%
% Then:  ic  = Pcc * icx
%        icx = pinv(Pcc) * ic
%        iv  = Pvv * ivx
%        ivx = pinv(Pvv) * iv
%        ix  = Pxx * ixx
%        ixx = pinv(Pxx) * ix
%
% Circuit-connected voltages, vxx: 
%        vxx = Pxx' * vx;
% 
% Circuit-connected resistances, rxx: 
%        rxx = diag(Pxx' * diag(rx) * Pxx);
%
% WARNING: pinv(Pxx) does not equal Pxx'
%
% ANOTHER WARNING: These transition maps assume the connected and unconnected
% current vectors are both in PER-TURN units. If not in per-turn units, one
% must first convert to per-turn using the turn #'s listed in, e.g., 
% tok_data_struct.ccnturn. 
%
% Lots of indices for indexing into specific vectors are also defined. See
% code comments.
%
% Josiah Wai 3/22

function circ = sparc_circ(tok_data_struct)

struct_to_ws(tok_data_struct);

% ================
% COIL GROUPINGS
% ================
% CS1 coils in series, VS coils antiseries, everything else independent
cccirc = [1 1 2 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 -19];

ccnames = cellstr(tok_data_struct.ccnames);       % coil names
cxnames = {'CS1U', 'CS1L', ccnames{5:20} 'VS1'};  % circuit names
ncx = length(cxnames);                            % num circuits

% The iy struct will hold indices of each coil/circuit so that they can be
% referred to easily by name (e.g, iy.CS1IU := 1)
for i = 1:nc
  iy.coils.(ccnames{i}) = i;
end

for i = 1:ncx
  iy.coil_circuits.(cxnames{i}) = i;
end


% Current fraction vectors
ccnturn = tok_data_struct.ccnturn';
ccfrac = zeros(size(ccnturn));
for icirc = 1:ncx
  icoils = find(abs(cccirc) == icirc);
  sgn = sign(cccirc(icoils));
  ccfrac(icoils) = sgn .* ccnturn(icoils) / sum(ccnturn(icoils));  
end


% ================
% VESSEL GROUPINGS
% ================

% vvgroup: group circuit elements according to the names defined in
% tok_data_struct.vvnames

vvdata = tok_data_struct.vvdata;  
vvnames = tok_data_struct.vvnames; 

vvnames = cellstr(vvnames);           % vessel names
vxnames = unique(vvnames, 'stable');  % vessel circuit/group names

igroup = 0;
for i = 1:length(vxnames)
  igroup = igroup + 1;
  j = find(strcmp(vvnames, vxnames{i}));
  vvgroup(j) = igroup;
end

vvcirc = 1:max(vvgroup);  % all vessel elements connected in parallel

% ================
% TRANSITION MAPS
% ================

% Pcc implements series connections for coil circuits (see header) 
Pcc = zeros(length(cccirc), max(cccirc));  
for i = 1:max(cccirc)
  Pcc(cccirc == i,i) = 1;
  Pcc(cccirc == -i,i) = -1;  
end

% Pvv implements parallel connections for vessel circuits (see header)
% For parallel connections, the current fraction is determined by vessel 
% resistances which follow 1/Rtot = 1/r1 + 1/r2 + 1/r3 ...
% Therefore current fraction in element 1 is: I1/Itot = (V/r1)/(V/R) = Rtot/r1
resv = tok_data_struct.resv;
Pvv = zeros(length(vvcirc), max(vvcirc)); 
for ii=1:max(vvcirc)
  idx1=find(vvgroup==ii);
  sum_rinv = sum(1./resv(idx1));
  vvfrac(idx1) = 1./resv(idx1)/sum_rinv;  
  Pvv(idx1,ii)=vvfrac(idx1);  
end

% combine circuit, vessel, plasma
Pxx = blkdiag(Pcc, Pvv, 1);

Pss = blkdiag(Pcc, Pvv);

% =======================================
% Define/copy various numbers and indices
% =======================================

nc = tok_data_struct.nc;  % num coils
nv = tok_data_struct.nv;  % num vessel elements
ns = nc + nv;             % num stabilizing elements
np = 1;                   % 1 plasma current
nx = nc + nv + np;        % total num of conducting elements

ncx = max(cccirc);        % num coil circuits
nvx = max(vvgroup);       % num vessel circuits
nsx = ncx + nvx;          % num stabilizing circuits
nxx = ncx + nvx + np;     % total num of circuits

iic = 1:nc;               % index of coils
iiv = (nc+1):(nc+nv);     % index of vessel elements
iis = 1:ns;               % index of stabilizing elements
iip = ns+1;               % index of plasma circuit
iix = 1:nx;               % index of all conducting elements

iicx = 1:ncx;             % index of coil circuits
iivx = (ncx+1):(ncx+nvx); % index of vessel circuits
iisx = 1:nsx;             % index of stabilizing circuits
iipx = nxx;                % index of plasma circuit
iixx = 1:nxx;              % index of all circuits


circ = variables2struct(ccnames, cxnames, vvnames, vxnames, iy, cccirc, ...
  ccfrac, vvcirc, vvgroup, Pcc, Pvv, Pss, Pxx, nc, nv, np, ns, nx, ncx, ...
  nvx, nsx, nxx, iic, iiv, iip, iis, iix, iicx, iivx, iisx, iipx, iixx);




















