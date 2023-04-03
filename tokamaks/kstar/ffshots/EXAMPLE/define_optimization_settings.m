% define various settings for the optimization including:
% 
%   timebase for optimization
%   which targets to explicitly control
%   voltage and current constraints

function settings = define_optimization_settings()

% time base for optimization
s.t0 = 0.4;
s.tf = 3;
s.Nlook = 30;
s.t = linspace(s.t0, s.tf, s.Nlook)';
s.dt = mean(diff(s.t));

% number of Grad-Shafranov iterations to perform
s.niter = 5;  

% number of vessel modes to use. 
% vessel modes are computed using a balanced realization on the vessel
% currents, see mpc_config.m
s.nvessmodes = 40;


% fds2control are the variables that are be explicitly controlled by
% the optimization algorithm. Also see measure_y and build_cmat. 
%
% ic                  - current in coils
% diff_psicp_psix     - flux error at control points vs x-point
% diff_psicp_psitouch - flux error at control points vs touch point
% psibry              - flux at boundary defining point
% psix_r              - flux derivative wrt r at target x-point
% psix_z              - flux derivative wrt z at target x-point

s.fds2control = {'ic', 'diff_psicp_psix', 'diff_psicp_psitouch', ...
  'psibry', 'psix_r', 'psix_z'}';


% all coils are active
s.active_coils = 1:13;  

% power supply voltage limits
s.enforce_voltage_limits = 0;
s.vmax = [1000 1000 500 500 1000 1000 1000 500 500 1000 1000 500 1000]';
s.vmin = -s.vmax;


% power supply current limits
s.enforce_current_limits = 0;
s.ic_max = [7.5 8.5 9 9 7 3.5 3.5 9 9 7 3.5 8 3]' * 1e3;
s.ic_min = -s.ic_max;

settings = s;



























