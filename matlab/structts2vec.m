function vec = structts2vec(s, fdnames, time)

x = cell(length(fdnames),1);
for i = 1:length(fdnames)  
  fd = fdnames{i};
  x{i} = interp1(s.(fd).Time, s.(fd).Data, time(:), 'linear', 'extrap'); 
end

x = horzcat(x{:})';

vec = x(:);









