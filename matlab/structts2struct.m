function snew = structts2struct(s, fdnames, time)

snew = struct;
for i = 1:length(fdnames)  
  fd = fdnames{i};
  if isstruct(s.(fd))
    snew.(fd) = interp1(s.(fd).Time, s.(fd).Data, time(:), 'linear', 'extrap');  
  end
end









