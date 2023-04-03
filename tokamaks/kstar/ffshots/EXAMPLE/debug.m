i = 15;
psi = reshape(y.psizr.Data(i,:), 65, 65);
ref = structts2struct(shapes, fields(shapes), t(i));
psicp = bicubicHermite(tok.rg, tok.zg, psi, ref.rb, ref.zb)';
psitouch = bicubicHermite(tok.rg, tok.zg, psi, ref.rtouch, ref.ztouch)';


x1 = psicp - psitouch;
x2 = y.diff_psicp_psitouch.Data(i,:)';


xk = [x.ic.Data(i,:) x.iv.Data(i,:)]';
C = config.cmats{i};
yk = yks{i} + C*xk;
x3 = yk(cv.iy.diff_psicp_psitouch);






