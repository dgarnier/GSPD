% Define time-dependent weights for each of the settings.fds2control (and
% also the actuator voltage). 
%
% weights.wts is the weight on the value of the parameter. weights.dwts is
% the weight on the derivative of the value. Parameters are not normalized
% so the weights can span orders of magnitude to find a good solution. 
%  
% To view a summary plot of the weights, set opts.plotit=1

function weights = define_optimization_weights(targs, settings, opts)

t = settings.t(:);
N = settings.Nlook;
ncp = size(targs.diff_psicp_psix.Data,2);   % number of shape control points
nu = length(settings.active_coils);    % number of actuators (power supplies)

% initialize all wts and dwts to zero
for dum = settings.fds2control(:)'
  fd = dum{:};
  wts.(fd).Time = t(:);
  dwts.(fd).Time = t(:);

  ny = size(targs.(fd).Data, 2);
  wts.(fd).Data  = zeros(N,ny);    
  dwts.(fd).Data = zeros(N,ny);    
end



% Populate the weights:

% weights on the shaping
wts.psibry.Data(:)           = 5e6;
wts.diff_psicp_psitouch.Data = sigmoidn(t, 0.15, 0.2, 1, 0) * ones(1,ncp) * 5e6;
wts.diff_psicp_psix.Data     = sigmoidn(t, 0.15, 0.2, 0, 1) * ones(1,ncp) * 5e6;
wts.psix_r.Data(:)       = sigmoidn(t, 0.15, 0.25, 0, 1) * 3e7;
wts.psix_z.Data(:)       = sigmoidn(t, 0.15, 0.25, 0, 1) * 3e7;


% weight the outer boundary point even higher
wts.diff_psicp_psitouch.Data(:,1) = wts.diff_psicp_psitouch.Data(:,1) * 30;
wts.diff_psicp_psix.Data(:,1) = wts.diff_psicp_psix.Data(:,1) * 30;


s = interp1([0 0.4 1], [5 1 1], t);
wts.diff_psicp_psix.Data = diag(s) * wts.diff_psicp_psix.Data;

% weights on control effort
wts.v.Time = t(:);
dwts.v.Time = t(:);
wts.v.Data = zeros(N,nu);
dwts.v.Data = ones(N,1) * [1 1 1 1 0.2 1 1 1];


weights = variables2struct(wts, dwts);


if opts.plotit
  fds = [settings.fds2control; 'v'];

  h = plot_structts(wts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.wts', 'fontsize', 14); drawnow

  h = plot_structts(dwts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.dwts', 'fontsize', 14); drawnow
end












