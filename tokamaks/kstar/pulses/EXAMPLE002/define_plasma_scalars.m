%  inputs: opts struct
%          if opts.plotlevel >= 2, makes a plot of the plasma scalars
% 
%  outputs: the plasma_scalars struct that has fields:
%     ip, li, wmhd, Rp each with subfields (Time, Data, Units). 
%     Units is only used for plotting, i.e. don't change units 
%     from A to kA and expect the code to intelligently adapt. 
%
function plasma_scalars   = define_plasma_scalars(opts)


% ip, plasma current in Amps
s.ip.Time = [0 0.65 1.3 10];
s.ip.Data = -[0 3.2 5.9 5.9] * 1e5;
s.ip.Units = 'A';

% li, internal inductance
s.li.Time = [0 0.6 2.2 2.75 3.5 10];
s.li.Data = [0.35 0.5 1.25 0.95 0.75 0.75];
s.li.Units = '';

% wmhd, stored thermal energy
s.wmhd.Time = [0 2.7 10];
s.wmhd.Data = [0 3.25 3] * 1e5;
s.wmhd.Units = 'J';

% Rp, plasma resistance in Ohms
s.Rp.Time = [0 1 3.9 10.1]';
s.Rp.Data = [1.4 0.7 0.3 0.3]' * 1e-6 * 0.4; 
s.Rp.Units = 'Ohm';

% format data
s = check_structts_dims(s); 

plasma_scalars = s;

% plotting
if opts.plotlevel >= 2
  plot_structts(s, fields(s), 2);
  sgtitle('Plasma scalar targets'); drawnow
end


































