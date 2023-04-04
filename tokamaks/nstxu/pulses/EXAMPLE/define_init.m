% define the initial condition for the optimization
%
% inputs  - none
% outputs - init struct
%           init.ic = initial coil currents
%           init.iv = initial vessel currents
%           init.psibry = initial boundary flux (used for Ip evolution)
%           init.v      = starting voltages for active coils
%
% For this example, we will grab these values from the presaved equilibrium
% from shot 204660 at t=30ms.

function init = define_init()

eq = load('eq204660_030.mat').eq;

init.ic = eq.icx;
init.iv = eq.ivx;
init.psibry = eq.psibry;

% initial voltages, only used to weight derivative of voltages
% can usually set to zero without affecting solution.
init.v  = zeros(8,1);  