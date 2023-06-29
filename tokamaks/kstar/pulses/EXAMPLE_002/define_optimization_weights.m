% Define time-dependent weights for each of the settings.fds2control (and
% also the actuator voltage). 
%
% weights.wts is the weight on the value of the parameter. weights.dwts is
% the weight on the derivative of the value. Parameters are not normalized
% so the weights may span orders of magnitude. 
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
  d2wts.(fd).Time = t(:);

  ny = size(targs.(fd).Data, 2);
  wts.(fd).Data  = zeros(N,ny);    
  dwts.(fd).Data = zeros(N,ny); 
  d2wts.(fd).Data = zeros(N,ny); 
end
wts.v.Time = t(:);
dwts.v.Time = t(:);
d2wts.v.Time = t(:);

wts.v.Data = zeros(N,nu);
dwts.v.Data = zeros(N,nu);
d2wts.v.Data = zeros(N,nu); 


%% Populate the weights:

% weight on psibry (for surface voltage and driving Ip)
wts.psibry.Data(:) = 5e6;  

% wt on flux err vs touch-point starts on and turns off as plasma diverts
wts.diff_psicp_psitouch.Data = sigmoidn(t, 1.5, 2.3, 1, 0) * ones(1,ncp) * 5e6;

% wt on flux err vs x-point starts off and turns on as plasma diverts
wts.diff_psicp_psix.Data  = sigmoidn(t, 2, 2.3, 0, 1) * ones(1,ncp) * 5e6;


% weight on flux gradient turns on as plasma diverts
wts.psix_r.Data(:) = sigmoidn(t, 1.5, 2, 0, 1) * 1e8;
wts.psix_z.Data(:) = sigmoidn(t, 1.5, 2, 0, 1) * 1e8;


% all coils are free except 12 and 13 (IRCUL and IVCUL) which are weighted
% (these will have a target of zero current)
wts.ic.Data  = ones(N,1) * [0 0 0 0 0 0 0 0 0 0 0 1e1 1e1];


% weight on the 2nd derivative of the coil currents (to penalize non-smooth
% trajectories)
d2wts.ic.Data = ones(N,1) * [1 1 1 1 1 1 1 1 1 1 1 1 1] * 0.2;


% no weight on absolute voltage
wts.v.Data = zeros(N,nu);   

% weight on the change in voltage - as is, the weights are roughly
% proportional to 1/coil inductance, but not too sensitive to the weighting
dwts.v.Data = ones(N,1) * [2 4 27 13 2 1 1 27 13 2 1 30 30] * 0.1;  


weights = variables2struct(wts, dwts, d2wts);


if opts.plotlevel >= 2
  fds = [settings.fds2control; 'v'];

  h = plot_structts(wts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.wts', 'fontsize', 14); drawnow

  h = plot_structts(dwts, fds, 3, [], 'linewidth', 1.5);
  sgtitle('weights.dwts', 'fontsize', 14); drawnow
end












