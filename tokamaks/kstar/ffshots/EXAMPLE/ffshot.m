clear all; clc; close all

%% Initialization
opts = struct;
opts.plotit = 0;
opts.debug = 0;
opts.verbose = 0;

tok              = load('kstar_tok').tok;
shapes           = define_shapes(opts);
plasma_scalars   = define_plasma_scalars(opts);
init             = define_init;
settings         = define_optimization_settings;
targs            = define_optimization_targets(shapes, tok, settings, opts);
weights          = define_optimization_weights(targs, settings, opts);



load('/Users/jwai/Research/GSpulse/tokamaks/kstar/externaldata/icdata.mat')
load('/Users/jwai/Research/GSpulse/tokamaks/kstar/externaldata/circ.mat')
x = resample(icdata, settings.t);
x.Data = squeeze(x.Data)';
y = x;
x.Data = x.Data * pinv(circ.Pcc)';
x.Data(:,[12 13]) = 0;
targs.ic = x;

%% Solve Grad-Shafanov + circuit dynamics
opts.plotit = 1;
soln = GSpulse(tok, shapes, plasma_scalars, init, settings, ...
  targs, weights, opts);



%% Plot results
summary_soln_plot(settings.t, shapes, soln.eqs, tok);

























