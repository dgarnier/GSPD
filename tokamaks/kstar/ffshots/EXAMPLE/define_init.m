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
% from shot 23436.

function init = define_init()


tok = load('kstar_tok.mat').tok;
efits = load('efits23436.mat').efits;


% init.ic = -[4.4 4.5 7 8.7 4 -0.4 -0.95 7 8.7 4 -0.4 2 0.17]' * 1e3;
init.ic = -[4.4 4.5 7 8.7 4 -0.4 -0.95 7 8.7 4 -0.4 0 0]' * 1e3;
% init.iv = zeros(tok.nv,1);
init.iv = load('init').init.iv;

init.v = zeros(tok.nc, 1);

% init.psibry = efits(1).psibry;
init.psibry = -2.7;



































