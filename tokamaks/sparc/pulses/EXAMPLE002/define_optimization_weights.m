% Define time-dependent weights for each of the settings.fds2control (and
% also the actuator voltage). 
%
% weights.wts is the weight on the value of the parameter. weights.dwts is
% the weight on the derivative of the value. Parameters are not normalized
% so the weights can span orders of magnitude to find a good solution. 
%  
% To view a summary plot of the weights, set opts.plotlevel >= 2

function weights = define_optimization_weights(targs, settings, opts)

% read some parameters
t = settings.t(:);
N = settings.N;
ncp = size(targs.diff_psicp_psix.Data,2); % number of shape control points
nu = length(settings.active_coils);       % number of actuators (power supplies)


% initialize all wts and dwts to zero
for dum = settings.fds2control(:)'
  fd = dum{:};
  wts.(fd).Time = t(:);
  dwts.(fd).Time = t(:);

  ny = size(targs.(fd).Data, 2);
  wts.(fd).Data  = zeros(N,ny);    
  dwts.(fd).Data = zeros(N,ny);    
end
wts.v.Time = t(:);
dwts.v.Time = t(:);
wts.v.Data = zeros(N,nu);
dwts.v.Data = zeros(N,nu);

%% Populate the weights:

% weight on psibry (for surface voltage and driving Ip)
wts.psibry.Data(:) = 5e6;  


% wt on flux err vs touch-point starts on and turns off as plasma diverts
wts.diff_psicp_psitouch.Data = sigmoidn(t, 2.5, 3, 1, 0) * ones(1,ncp) * 5e6;


% wt on flux err vs x-point starts off and turns on as plasma diverts
wts.diff_psicp_psix.Data  = sigmoidn(t, 2.5, 3, 0, 1) * ones(1,ncp) * 5e6;


% weight on flux gradient turns on as plasma diverts
wts.psix_r.Data(:) = sigmoidn(t, 2.5, 3, 0, 1) * 3e7;
wts.psix_z.Data(:) = sigmoidn(t, 2.5, 3, 0, 1) * 3e7;


% weight the outer boundary point even higher
% wts.diff_psicp_psitouch.Data(:,1) = wts.diff_psicp_psitouch.Data(:,1) * 30;
% wts.diff_psicp_psix.Data(:,1)     = wts.diff_psicp_psix.Data(:,1) * 30;


% all coils free except VSC
% wts.ic.Data  = ones(N,1) * [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1e1];
wts.ic.Data  = ones(N,1) * [ones(1,18)*1e-3 1];


% no weight on absolute voltage
wts.v.Data = zeros(N,nu);   

% weight on the change in voltage - as is, the weights are roughly
% proportional to 1/coil inductance
dwts.v.Data = ones(N,1) * [1 1 3 3 3 3 4 4 1 1 1 1 0.5 0.5 20 20 20 20 100] * 0.1;


weights = variables2struct(wts, dwts);


if opts.plotlevel >= 2
  fds = [settings.fds2control; 'v'];

  h = plot_structts(wts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.wts', 'fontsize', 14); drawnow

  h = plot_structts(dwts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.dwts', 'fontsize', 14); drawnow
end












