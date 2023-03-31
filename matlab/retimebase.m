function X = retimebase(X, t)


% interpolate onto new timebase
fds = fields(X);
for i = 1:length(fds)
  fd = fds{i};

  if isstruct(X.(fd)) && isfield(X.(fd), 'Data')
    X.(fd).Data = interp1hold(X.(fd).Time, X.(fd).Data, t);
    X.(fd).Time = t;
  end
end
