% Creates a ui figure for plotting equilibrium solutions

function fig = summary_soln_plot(times, shapes, eqs, tok)

% create figure
fig = figure;
fig.Position = [754   322   425   611];
ax = axes(fig);
ax.Box = 'on';
ax.Position = [0.13 0.2 0.775 0.7];


% slider button
s = uicontrol(fig, 'style', 'slider');
s.Units = 'normalized';
s.Position = [0.15 0.05 0.6 0.05];
s.Min = times(1);
s.Max = times(end);
s.Value = times(1);
s.Callback = {@sliderCallback, times, shapes, eqs, tok};

% text edit button
e = uicontrol(fig, 'style', 'edit');
e.Units = 'normalized';
e.Position = [0.8 0.07 0.15 0.05];
e.Callback = {@editCallback, times, shapes, eqs, tok};

plot_shape(s.Value, times, shapes, eqs, tok)


% slider callback
function sliderCallback(src, event, times, shapes, eqs, tok)
  t = src.Value;
  plot_shape(t, times, shapes, eqs, tok)
end


% text edit callback
function editCallback(src, event, times, shapes, eqs, tok)
  t = str2double(src.String);
  plot_shape(t, times, shapes, eqs, tok)
end


% plot shape targets
function plot_shape(t, times, shapes, eqs, tok)

  [~,i] = min(abs(t-times));

  ref = structts2struct(shapes, fields(shapes), t);

  cla
  hold on
  % eq = find_bry(eqs{i}.psizr, tok, 0);
  plot_eq(eqs{i}, tok, 'r', 'linewidth', 1.5)
  scatter(ref.rb, ref.zb, 20, 'b', 'filled')
  plot(ref.rx, ref.zx, 'xb', 'linewidth', 4, 'markersize', 14)
  scatter(ref.rtouch, ref.ztouch, 100, 'db', 'filled')
  text(-0.25, -0.1, 'Drag slider to view shape targets.', ...
    'units', 'normalized', 'fontsize', 11)
  str = sprintf('eq %d: time=%.3f', i, t);
  title(str, 'fontsize', 14)
  drawnow

end


end





































