% define the initial condition for the optimization
%
% inputs  - none
% outputs - init struct
%           init.ic = initial coil currents
%           init.iv = initial vessel currents          
%           init.v  = initial voltages in the active coils, only used for
%                     weighting the derivative of voltages and can often
%                     be set to zero. 
%

function init = define_init()

% these numbers for coil currents are from shot 23436 at t=0.4
init.ic = -[4.4 4.5 7 8.7 4 -0.4 -0.95 7 8.7 4 -0.4 0 0]' * 1e3;
init.iv = load('init').init.iv;
init.v = zeros(13, 1);






























