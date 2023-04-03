function mpcsoln = mpc_update_psiapp(pcurrt, config, tok, shapes, ...
  plasma_scalars, init, settings, targs, weights, opts)

% read parameters
fds2control = settings.fds2control;
Nlook = settings.Nlook;
t = settings.t;
dt = settings.dt;
cv = config.cv;
nu = config.nu;
nx = config.nx;
E = config.E;
F = config.F;
Fw = config.Fw;
Chat = config.Chat;
wts = weights.wts;
dwts = weights.dwts;


% target psibry that is consistent with targ.ip
psipla = tok.mpp * pcurrt;
% psipla0 = reshape(psipla(:,1), tok.nz, tok.nr);
% psiapp0 = tok.mpc*init.ic + tok.mpv*init.iv;
% psiapp0 = reshape(psiapp0, tok.nz, tok.nr);
% psi0 = psiapp0 + psipla0;
% eq0 = find_bry(psi0, tok, 0);
% init.psibry = eq0.psibry;
init.psibry = -2.25;

targs.psibry = psibry_dynamics(tok, settings, shapes, plasma_scalars,...
  init, psipla);


% measure y
yks = measure_ys(psipla, fds2control, shapes, tok);
ykhat = vertcat(yks{:});

% target y and error 
rhat = structts2vec(targs, fds2control, t);
dytarghat = rhat - ykhat;
ny = length(yks{1});


% initial state
uprev = init.v;
xk = [init.ic; init.iv];
xk = config.bal.T * xk;   
x0 = zeros(size(xk));
dxk = xk;
x0hat = repmat(x0, Nlook, 1);


% plasma-coupling term
w = plasma_coupling(settings.dt, tok, pcurrt);
w = config.bal.T * w;
[~,wd] = c2d(config.Ar, w, dt);
wd = wd(:);


% weights
q = structts2vec(wts, fds2control, t);
r = structts2vec(wts, {'v'}, t);
dr = structts2vec(dwts, {'v'}, t);


% filter nans
i = isnan(dytarghat);
dytarghat(i) = 0;
q(i) = 0;
Chat(i,:) = 0;

% Note: ehat = M*duhat + d
M = -Chat*F;
d = dytarghat - Chat * (E*dxk + Fw*wd);

% sparse format weight matrices
Q = spdiags(q, 0, length(q), length(q));
R = spdiags(r, 0, length(r), length(r));
dR = spdiags(dr, 0, length(dr), length(dr));

% J1
H1 = M'*Q*M;
f1 = M'*Q*d;

% J3
H3 = R;
f3 = 0;

% J4
m = tridiag(-1, 1, 0, Nlook);
Su = kron(m, eye(nu));
Su = sparse(Su);

H4 = Su' * dR * Su;
f4 = -Su' * dR(:,1:nu) * uprev;


% J total
f = f1 + f3 + f4;
H = H1 + H3 + H4;
H = (H+H')/2;

% no equality constraints
Aeq = [];
beq = [];

% Inequality constraints

% voltage limits
if settings.enforce_voltage_limits
  ub = repmat(settings.vmax, Nlook, 1);
  lb = repmat(settings.vmin, Nlook, 1);
else
  ub = [];
  lb = [];
end

% current limits
if settings.enforce_current_limits
  ymin = -inf(ny,1);
  ymax = inf(ny,1);
  ymin(cv.iy.ic) = settings.ic_min;
  ymax(cv.iy.ic) = settings.ic_max;
  yminhat = repmat(ymin, Nlook, 1);
  ymaxhat = repmat(ymax, Nlook, 1);
  Aineq = [-M; M];
  bineq = [ymaxhat + d - rhat; -yminhat - d + rhat];
  i = isinf(bineq);
  Aineq(i,:) = [];
  bineq(i,:) = [];
else
  Aineq = [];
  bineq = [];
end


% solve quadratic program
duhat0 = -H\f;
qpopts = optimoptions('quadprog', 'algorithm', 'active-set', 'Display', 'off');
duhat = quadprog(H,f,Aineq, bineq, Aeq, beq, lb, ub, duhat0, qpopts);


% extract predictions
ehat = M*duhat + d;
dxhat = E*dxk + F*duhat;
dyhat = dytarghat - ehat;
yhat = dyhat + ykhat;
xhat = dxhat + x0hat;
    
y = vec2slstructts(yhat, fds2control, cv.iy, t);
x = vec2slstructts(xhat, fields(cv.ix), cv.ix, t);

y = copyfields(y,x,[],0);

psiapp = [tok.mpc tok.mpv] * config.bal.Ti * [x.ic.Data'; x.iv.Data'];
% psiapp = tok.mpc * x.ic.Data' + tok.mpv * x.iv.Data';

psizr = psiapp + psipla;

y.psizr.Time = t;
y.psiapp.Data = t;
y.psipla.Data = t;

y.psizr.Data = psizr';
y.psiapp.Data = psiapp';
y.psipla.Data = psiapp';

mpcsoln = y;

% plot stuff, debugging
if opts.plotit
  h = plot_structts(y, fds2control, 2);
  h = plot_structts(targs, fds2control, 2, h, '--r');
  drawnow

 
  ic = targs.ic.Data';
  xr = vec2slstructts(ic(:), tok.ccnames, cv.ix, t);
  h = plot_structts(x, tok.ccnames, 2);
  h = plot_structts(xr, tok.ccnames, 2, h, '--r');
  drawnow

end


if 0
  i = 50;
  figure
  psi = reshape(psizr(:,i), tok.nz, tok.nr);
  eq = find_bry(psi, tok, 0);
  plot_eq(eq, tok, 'r')
  title(settings.t(i))  
end





































