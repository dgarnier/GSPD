% builds dynamic model and some MPC stuff that can be precomputed

function config = mpc_config(tok, shapes, targs, settings)


% cv will hold indices of the outputs (y), states (x), and inputs (u). 
% Useful for organizing data.
cv = struct;
cv.ynames = settings.fds2control;
idx = 0;
for k = 1:length(cv.ynames)    
  varname = cv.ynames{k};
  n = size(targs.(varname).Data, 2);  
  idx = idx(end)+1:idx(end)+n;  
  cv.iy.(varname) = idx;
end 
cv.xnames = {'ic', 'iv'}';
cv.ix.ic = 1:tok.nc;
cv.ix.iv = tok.nc + (1:settings.nvessmodes);
for i = 1:tok.nc
  cv.ix.(tok.ccnames{i}) = i;
  cv.iu.(tok.ccnames{i}) = i;
end
cv.unames = {'v'};
nu = length(settings.active_coils);
cv.iu.v = 1:nu;


% build the dynamics model A,B matrices
M = [tok.mcc tok.mcv; tok.mcv' tok.mvv];
M = (M + M') /  2;
R = diag([tok.resc; tok.resv]);                   
Minv = inv(M);
A = -Minv * R;
B = Minv(:,settings.active_coils);


% build the output C matrices
dpsizrdx = [tok.mpc tok.mpv];
cmats = output_model(dpsizrdx, tok, shapes, settings);


% perform a balanced realization on the vessel currents
% step1: compute balancing transformation matrices
ivess = tok.nc + (1:tok.nv);
Avess = A(ivess,ivess);
Bvess = B(ivess,:);
Cvess = eye(tok.nv, tok.nv);
Pvess = ss(Avess,Bvess,Cvess,0);
[bsys,g,T,Ti] = balreal(Pvess);
T = blkdiag(eye(tok.nc), T);
Ti = inv(T);
iuse = 1:(tok.nc + settings.nvessmodes);
T = T(iuse,:);
Ti = Ti(:,iuse);

bal.T = T;
bal.Ti = Ti;
bal.info = ['Balanced realization transforming state vector x=[ic;iv] to ' ...
  'xb=[ic; ivb]. The transformations are xb=T*x and x=Ti*xb.'];


% step2: reduce models
for i = 1:length(cmats)
  cmats{i} = cmats{i} * Ti;
end
Chat = blkdiag(cmats{:});
Ar = T*A*Ti;
Br = T*B;
[nx, nu] = size(Br);


% discretize model
[Ad, Bd] = c2d(Ar,Br,settings.dt);


% Prediction model used in MPC. See published paper for definitions.
Nlook = settings.Nlook;
nw = nx;
E = [];
F  = [];
Fw = [];
Apow  = eye(nx);
F_row = zeros(nx, Nlook*nu);
Fw_row = zeros(nx, Nlook*nw);

for i = 1:Nlook

  idx = (nu*(i-1)+1):(nu*i);
  F_row = Ad * F_row;
  F_row(:,idx) = Bd;
  F = [F; F_row];

  idx = (nw*(i-1)+1):(nw*i);
  Fw_row = Ad * Fw_row;
  Fw_row(:,idx) = eye(nx);
  Fw = [Fw; Fw_row];
  
  Apow = Ad * Apow;
  E = [E; Apow];
end


% save prediction model
config = variables2struct(cv,nx,nu,nw,A,B,Ar,Br,Ad,Bd,Minv,E,F,Fw,Chat,...
  cmats,bal);








































































