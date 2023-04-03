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
init.ic = -[4.4 4.5 7 8.7 4 -0.4 -0.95 7 8.7 4 -0.4 0 0]' * 1e3;
init.iv = load('init').init.iv;
init.v = zeros(tok.nc, 1);
init.psibry = -2.7;



% cccirc = [1 2 3 4 5 6 7 1 2 8 9 10 11 7 13 12 -13 12];
% Pcc = cccirc_to_Pcc(cccirc);
% eq = load('init').init; 
% init.ic = pinv(Pcc) * eq.ic;
% init.iv = eq.iv;
% init.psibry = eq.psibry;
% init.v = zeros(13,1);































