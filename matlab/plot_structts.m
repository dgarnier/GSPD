
function fighandle = plot_structts(s, fdnames, nrows, fighandle, varargin)

N = length(fdnames);
if ~exist('nrows','var') || isempty(nrows), nrows = min(3,N); end
if ~exist('fighandle','var') || isempty(fighandle), fighandle = figure; end

ncols = ceil(N/nrows);

figure(fighandle);

for i = 1:N
  fd = fdnames{i};
  subplot(nrows, ncols, i)
  hold on
  grid on
  plot( s.(fd).Time, s.(fd).Data, varargin{:});
  title(fd, 'fontsize', 16, 'fontweight', 'bold')
  if isfield(s.(fd), 'Units')
    ylabel(s.(fd).Units, 'fontsize', 12)
  end
  xlabel('time')
end

  






























