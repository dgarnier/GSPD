
function s = vec2slstructts(vec, fdnames, iy, time)

x = reshape(vec, [], length(time))';

for i = 1:length(fdnames)  
  fd = fdnames{i};

  s.(fd).Time = time(:);
  s.(fd).Data = x(:, iy.(fd));
end

  