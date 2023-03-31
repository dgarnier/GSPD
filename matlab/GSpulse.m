function soln = GSpulse(tok, shapes, plasma_scalars, init, settings, ...
  targs, weights, opts)


% put everthing on same timebase
t = settings.t;
shapes         = retimebase(shapes, t);
targs          = retimebase(targs, t);
plasma_scalars = retimebase(plasma_scalars, t);
weights.wts    = retimebase(weights.wts, t);
weights.dwts   = retimebase(weights.dwts, t);


% dynamics model and any mpc stuff that can be precomputed 
config = mpc_config(tok, shapes, targs, settings);  


% initialize with estimate of plasma current distribution
pcurrt = initialize_pcurrt(tok, shapes, plasma_scalars);
psipla = tok.mpp * pcurrt;


% perform iterations between dynamics optimization and Grad Shafranov 
args = {config, tok, shapes, plasma_scalars, init, settings, targs, ...
  weights, opts};

for iter = 1:settings.niter

  fprintf('\n\n\nGrad-Shafranov iteration %d of %d\n\n\n', iter, settings.niter)

  % dynamics optimization
  mpcsoln = mpc_update_psiapp(pcurrt, args{:});
  
  % Grad-Shafranov iteration
  [eqs, pcurrt] = gs_update_psipla(mpcsoln, tok, ...
    plasma_scalars, settings);

end


soln.eqs = eqs;
soln.mpcsoln = mpcsoln;



















