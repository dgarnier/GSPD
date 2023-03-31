
function vec = struct2vec(s, fdnames)

x = cell(length(fdnames),1);
for i = 1:length(fdnames)    
  x{i} = s.(fdnames{i});
end
vec = vertcat(x{:});



