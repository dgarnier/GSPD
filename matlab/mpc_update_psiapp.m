% Sets up and solves an MPC-like quadratic program for creating
% Grad-Shafranov equilibria according to the target shape evolution. 
% 
% Inputs: see GSpulse.m and pulse.m in the EXAMPLE folder
% Outputs: mpcsoln - struct with timeseries data on the currents and flux
%                    errors and flux evolution          

function mpcsoln = mpc_update_psiapp(iter, pcurrt, config, tok, shapes, ...
  plasma_scalars, init, settings, targs, weights, opts)


% read parameters
fds2control = settings.fds2control;
N = settings.N;
t = settings.t;
dt = settings.dt;
wts = weights.wts;
dwts = weights.dwts;
cv = config.cv;
nu = config.nu;
nx = config.nx;
E = config.E;
F = config.F;
Fw = config.Fw;
Cmats = config.Cmats;
Chat = blkdiag(Cmats{:});


% map plasma current -> plasma flux
psipla = tok.mpp * pcurrt;


% find the targ.psibry that is consistent with targ.ip
psiapp0 = tok.mpc*init.ic + tok.mpv*init.iv;
psi0 = psiapp0 + psipla(:,1);
psi0 = reshape(psi0, tok.nz, tok.nr);
ref = structts2struct(shapes, {'rb','zb'}, t(1));
psibry0 = mean(bicubicHermite(tok.rg, tok.zg, psi0, ref.rb, ref.zb));
targs.psibry = psibry_dynamics(tok, settings, shapes, plasma_scalars,...
  psibry0, psipla);


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
xk = config.bal.Tx * xk;   
x0 = zeros(size(xk));
dxk = xk;
x0hat = repmat(x0, N, 1);


% plasma-coupling term
w = plasma_coupling(settings.dt, tok, pcurrt);
w = config.bal.Tx * w;
[~,wd] = c2d(config.Ar, w, dt);
wd = wd(:);


% weights
q = structts2vec(wts, fds2control, t);
r = structts2vec(wts, {'v'}, t);
dr = structts2vec(dwts, {'v'}, t);

% filter out any nans
i = isnan(dytarghat);
dytarghat(i) = 0;
q(i) = 0;
Chat(i,:) = 0;


% form weight matrices, sparse format
Q = spdiags(q, 0, length(q), length(q));
R = spdiags(r, 0, length(r), length(r));
dR = spdiags(dr, 0, length(dr), length(dr));


% Prediction model and cost matrices (note: ehat = M*duhat + d)
% see accompanying paper for definitions and derivations
M = -Chat*F;
d = dytarghat - Chat * (E*dxk + Fw*wd);

% J1
H1 = M'*Q*M;
f1 = M'*Q*d;

% J3
H3 = R;
f3 = 0;

% J4
m = tridiag(-1, 1, 0, N);
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
  ub = repmat(settings.vmax, N, 1);
  lb = repmat(settings.vmin, N, 1);
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
  yminhat = repmat(ymin, N, 1);
  ymaxhat = repmat(ymax, N, 1);
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
dxhat = E*dxk + F*duhat + Fw*wd;
dyhat = dytarghat - ehat;
yhat = dyhat + ykhat;
xhat = dxhat + x0hat;
    
y = vec2structts(yhat, fds2control, cv.iy, t);
x = vec2structts(xhat, fields(cv.ix), cv.ix, t);
y = copyfields(y,x,[],0);

psiapp = [tok.mpc tok.mpv] * config.bal.Txi * [x.ic.Data'; x.ivb.Data'];
psizr = psiapp + psipla;

y.psizr.Time = t;
y.psiapp.Data = t;
y.psipla.Data = t;

y.psizr.Data = psizr';
y.psiapp.Data = psiapp';
y.psipla.Data = psiapp';

y.iv = y.ivb;
y.iv.Data = y.ivb.Data * config.bal.Tvi';

mpcsoln = y;

% plot timetraces
if (opts.plotlevel >= 2) || (opts.plotlevel >= 1 && iter==settings.niter)
  h = plot_structts(y, fds2control, 2);
  plot_structts(targs, fds2control, 2, h, '--r');
  legend('Actual', 'Target', 'fontsize', 16)
  drawnow
 
  ic = targs.ic.Data';
  xr = vec2structts(ic(:), tok.ccnames, cv.ix, t);
  h = plot_structts(x, tok.ccnames, 2);
  plot_structts(xr, tok.ccnames, 2, h, '--r');
  legend('Actual', 'Target', 'fontsize', 16)
  drawnow
end





































