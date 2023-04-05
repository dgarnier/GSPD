% Create the tok object from the external data sources (i.e. toksys/MEQ
% derived sources) for sparc. 

clear; clc; close all

saveit = 1;
savefn = [getenv('GSROOT') '/tokamaks/sparc/tok/sparc_tok.mat'];

% load geometry
tokdata = load('sparc_obj_6565.mat').tok_data_struct;
circ = sparc_circ(tokdata);
Pcc = circ.Pcc;

% mutual inductances and resistances
resc = diag(Pcc' * diag(tokdata.resc) * Pcc);
resv = tokdata.resv;
mcc = Pcc' * tokdata.mcc * Pcc;
mcv = Pcc' * tokdata.mcv;
mvv = tokdata.mvv;
mpc = tokdata.mpc * Pcc;
mpv = tokdata.mpv;
mpp = tokdata.mpp; 
ccnames = circ.cxnames(:);
nc = circ.ncx;
nv = circ.nv;
rg = tokdata.rg;
zg = tokdata.zg;
nr = tokdata.nr;
nz = tokdata.nz;
limdata = tokdata.limdata;

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
































