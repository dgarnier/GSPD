clear all; clc; close all

%% Initialization
opts = struct;
opts.plotit = 0;
opts.debug = 0;
opts.verbose = 0;

tok              = load('nstxu_tok').tok;
shapes           = define_shapes(opts);
plasma_scalars   = define_plasma_scalars(opts);
init             = define_init;
settings         = define_optimization_settings;
targs            = define_optimization_targets(shapes, tok, settings, opts);
weights          = define_optimization_weights(targs, settings, opts);


%% Solve Grad-Shafanov + circuit dynamics
soln = GSpulse(tok, shapes, plasma_scalars, init, settings, ...
  targs, weights, opts);



%% Plot results
summary_soln_plot(settings.t, shapes, soln.eqs, tok);

























