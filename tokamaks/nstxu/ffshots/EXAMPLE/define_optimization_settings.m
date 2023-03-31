% define various settings for the optimization including:
% 
%   timebase for optimization
%   which targets to explicitly control
%   voltage and current constraints

function settings = define_optimization_settings()

% time base for optimization
s.t0 = 0.05;
s.tf = 1;
s.Nlook = 50;
s.t = linspace(s.t0, s.tf, s.Nlook)';
s.dt = mean(diff(s.t));

% number of Grad-Shafranov iterations to perform
s.niter = 5;  

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


% only some of the coils are active on nstxu, corresponds to tok.ccnames
s.active_coils = [1 2 5 6 8 9 10 13];  

% power supply voltage limits
s.enforce_voltage_limits = 1;
s.vmax = [4048 1012 2024 2024 3036 2024 2024 1012]';
s.vmin = -s.vmax;


% power supply current limits
s.enforce_current_limits = 1;
s.ic_max = [20  15 inf  inf 15 8   inf   0   8 15 inf inf 15]' * 1e3;
s.ic_min = [-20 0 -inf -inf 0 -13 -inf -24 -13 0 -inf -inf 0]' * 1e3;


settings = s;



























