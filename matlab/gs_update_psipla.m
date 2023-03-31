function [eqs, pcurrtdata] = gs_update_psipla(...
  mpcsoln, tok, plasma_scalars, settings)


% read inputs
N = settings.Nlook;
psizr = mpcsoln.psizr.Data';
psiapp = mpcsoln.psiapp.Data';
nz = tok.nz;
nr = tok.nr;
mpp = tok.mpp;


% trace boundary
eqs = cell(N,1);
for i = 1:N
  fprintf('Tracing boundary: %d of %d ...\n', i, N);
  psizr_i = reshape(psizr(:,i), tok.nz, tok.nr);
  eqs{i} = find_bry(psizr_i, tok, 0);
end


if 0
  close all
  i = 11;
  plot_eq(eqs{i}, tok, 'r', 'linewidth', 1)
  scatter(shapes.rb.Data(i,:), shapes.zb.Data(i,:), 'k', 'filled')
  scatter(shapes.rx.Data(i), shapes.zx.Data(i), 100, 'b', 'filled')
end
  

pcurrtdata = zeros(nz*nr,N);
for i = 1:N

  eq = eqs{i};
  ip = plasma_scalars.ip.Data(i);
  wmhd = plasma_scalars.wmhd.Data(i);
  li = plasma_scalars.li.Data(i);
  
  psibry = eq.psibry;
  psimag = eq.psimag;
  
  [rgg, zgg] = meshgrid(tok.rg, tok.zg);
  dr = mean(diff(tok.rg));
  dz = mean(diff(tok.zg));
  dA = dr*dz;
  
  in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
  psinzr = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);
  psinzr(~in) = nan;
  
  mu0 = pi*4e-7;
  psin = linspace(0,1,tok.nr)';
  
  % will hold all the basis functions
  b = struct; 
  b.psin = psin;
  
  b.p1 = -4 * (-(psin-0.5).^2 + 0.25);                  % basis function for P'
  b.f1 = -polyval([0.54 -0.08 -1.46 1], psin) * 1e-6;   % 1st basis fun for FF'
  b.f2 = -ones(size(psin)) * 1e-6;                      % 2nd basis fun for FF'
  
  % pressure basis
  b.pres = cumtrapz(psin, b.p1);  
  b.pres = b.pres-b.pres(end);
  b.preszr = interp1(psin, b.pres, psinzr(:));
  in = ~isnan(b.preszr);
  b.preszr(~in) = 0;
  
  % pprime basis
  b.pprime = b.p1 * 2*pi/(psibry-psimag);
  b.pprimezr = interp1(psin, b.p1, psinzr(:)) * 2*pi/(psibry-psimag);
  b.pprimezr(~in) = 0;
  
  % ffprim basis
  b.ffprim = [b.f1 b.f2] * 2*pi/(psibry-psimag);
  b.ffprimzr = interp1(psin, [b.f1 b.f2], psinzr(:)) * 2*pi/(psibry-psimag);
  b.ffprimzr([~in ~in]) = 0;
  
  % set up the equations:  [Ip; wmhd; li] = H * [cp1; cf1; cf2]; 
  % then this can be solved for [cp1; cf1; cf2] the P' and FF' coefficients
  H = zeros(3);
  R = rgg(:);
  
  H(1,:) = [R'*b.pprimezr  1./(mu0 * R') * b.ffprimzr] * dA;
  H(2,1) = 3*pi*R'*b.preszr*dA;
  
  alpha = interp1([0.5 1.2], [-4 -10], li, 'linear', 'extrap');
  H(3,:) = [0 1 alpha];
  
  c = H \ [ip; wmhd; 0];
  
  jphi = [R.*b.pprimezr  1./(mu0 * R).*b.ffprimzr] * c;
  jphi = reshape(jphi, nz, nr);
  pcurrt = jphi * dA;
  
  
  eq.pcurrt = pcurrt;
  eq.psiapp = reshape(psiapp(:,i), nz, nr);
  eq.psipla = mpp * pcurrt(:);
  eq.psipla = reshape(eq.psipla, nz, nr);
  eq.psizr = eq.psiapp + eq.psipla;
  

  eqs{i} = eq;
  pcurrtdata(:,i) = pcurrt(:);
end



if 0
  figure
  contourf(tok.rg, tok.zg, pcurrt); colorbar
  set(gcf, 'Position', [731 342 408 624]);
end

if 0
  [Br, Bz, Bp] = psizr2b(eq.psizr, tok);
  volumezr = 2 * pi * (tok.rgg+dr/2) * dA;
  volumezr(~in) = 0;
  Vtot = sum(sum(volumezr));
  Bp2volavg = sum(sum((Br.^2 + Bz.^2).*volumezr)) / Vtot;
  Cl =  sum(sqrt(diff(eq.rbbbs).^2 + diff(eq.zbbbs).^2));
  Bp2bryavg = (mu0*ip/Cl)^2;
  li = Bp2volavg / Bp2bryavg;
end
























