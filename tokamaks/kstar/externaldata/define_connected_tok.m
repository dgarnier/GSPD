% Create the connected tok_data_struct object from the external data
% sources for kstar. 
clear; clc; close all

saveit = 1;
savefn = [getenv('GSROOT') '/tokamaks/kstar/tok/kstar_tok.mat'];

% load a bunch of shit
circ = load('circ').circ;
Pcc = circ.Pcc;
ktok = load('kstar_obj_mks_struct_2015_6565.mat').tok_data_struct;

resc = diag(Pcc' * diag(ktok.resc) * Pcc);
resv = ktok.resv;
ccnames = load('PS').PS.names;
mcc = Pcc' * ktok.mcc * Pcc;
mcv = Pcc' * ktok.mcv;
mvv = ktok.mvv;
mpc = ktok.mpc * Pcc;
mpv = ktok.mpv;
nc = length(ccnames);
nv = ktok.nv;
rg = ktok.rg;
zg = ktok.zg;
nr = ktok.nr;
nz = ktok.nz;
limdata = ktok.limdata';
mpp = ktok.mpp; 
% mpp = unwrap_mpp(mpp, ktok.nz, ktok.nr);
% mpp = (mpp+mpp')/2;


% descriptions
d.mcc = 'mutual inductances from coils to coils [Wb/A]';
d.mcv = 'mutual inductances from coils to vessels [Wb/A]';
d.mvv = 'mutual inductances from vessels to vessels [Wb/A]';
d.mpc = 'mutual inductances from plasma grid to coils [Wb/A]';
d.mpv = 'mutual inductances from plasma grid to vessels [Wb/A]';
d.mpp = 'see unwrapp_mpp.m, mutual inductances from plasma grid to plasma grid [Wb/A]';
d.resc = 'coil resistances [Ohms]';
d.resv = 'vessel resistances [Ohms]';
d.ccnames = 'coil names';
d.nc = 'number of coils';
d.nv = 'number of vessel elements';
d.rg = 'radial positions on grid [m]';
d.zg = 'vertical positions on grid [m]';
d.nr = 'grid size (radial)';
d.nz = 'grid size (vertical)';
d.limdata = '(z,r) of vertices defining limiter';
descriptions = d;

% save data
tok = variables2struct(mcc, mcv, mvv, mpc, mpv, mpp, resc, resv, ...
  ccnames, nc, nv, rg, zg, nr, nz, limdata); 

fds = sort(fields(tok));
tok = reorderstructure(tok, fds{:});
descriptions = reorderstructure(descriptions, fds{:}); 
tok.descriptions = descriptions; 

if saveit
  save(savefn, 'tok')
end
































