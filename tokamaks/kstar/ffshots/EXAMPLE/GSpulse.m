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
% init.v = -pinv(config.B) * config.A * [init.ic; init.iv];
% init.v = -tok.mcv*(inv(tok.mvv))*diag(tok.resv)*init.iv + tok.resc.*init.ic;

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
  [eqs1, eqs0, pcurrt1] = gs_update_psipla(mpcsoln, tok, shapes,...
    plasma_scalars, settings);

  a = 1;
  pcurrt = a*pcurrt1 + (1-a)*pcurrt;

  if iter == 2
    summary_soln_plot(settings.t, shapes, eqs0, tok);
    pause
  end

%   i = 2;
%   close all
%   eq = find_bry(eqs{i}.psizr, tok, 0);
%   plot_eq(eq, tok, 'r', 'linewidth', 1)
%   scatter(shapes.rb.Data(i,:), shapes.zb.Data(i,:), 'k', 'filled')
%   scatter(shapes.rx.Data(i), shapes.zx.Data(i), 100, 'b', 'filled')

end


soln.eqs = eqs1;
soln.mpcsoln = mpcsoln;



for i = 1:settings.Nlook
  disp(i)
  eqs1{i} = find_bry(eqs1{i}.psizr, tok, 0);
end
summary_soln_plot(settings.t, shapes, eqs1, tok);
















