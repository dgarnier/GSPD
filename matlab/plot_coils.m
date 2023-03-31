% Helper plotting function. 
%
% Makes a plot of coil values (such as coil currents, or coil voltages). 
% Plots them individually for each coil grouping (CS, PF, DV, VS)


function ax = plot_coils(t, y, tok_data_struct, ystring, varargin)

circ = sparc_circ(tok_data_struct);
iy = circ.iy.coil_circuits;

if size(y,1) ~= length(t)
  y = y'; 
end

ax(1) = subplot(221);
i = iy.CS1U:2:iy.CS3L;
plot(t, y(:,i), varargin{:})
legend(circ.cxnames{i}, 'location', 'best', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
ylabel(ystring, 'fontsize', 14)


ax(2) = subplot(222);
i = iy.PF1U:2:iy.PF4L;
plot(t, y(:,i), varargin{:})
legend(circ.cxnames{i}, 'location', 'best', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
ylabel(ystring, 'fontsize', 14)


ax(3) = subplot(223);
i = iy.DV1U:2:iy.DV2L;
plot(t, y(:,i), varargin{:})
legend(circ.cxnames{i}, 'location', 'best', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
ylabel(ystring, 'fontsize', 14)



ax(4) = subplot(224);
i = iy.VS1;
plot(t, y(:,i), varargin{:})
legend(circ.cxnames{i}, 'location', 'best', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
ylabel(ystring, 'fontsize', 14)
ylim([-1 1])

linkaxes(ax, 'x')

