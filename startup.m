% the startup file just adds various paths

GSROOT = fileparts(mfilename('fullpath'));
setenv('GSROOT', GSROOT);

addpath([GSROOT '/matlab'])


% select a tokamak
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath([GSROOT '/tokamaks']))
warning('on','MATLAB:rmpath:DirNotFound')
d = dir([GSROOT '/tokamaks']);
d(1:2) = [];
d(~[d(:).isdir]) = [];
n = length(d);
fprintf('Choose a tokamak [default, 1]: \n')
for i = 1:n
  fprintf('[%d] %s \n', i, d(i).name)
end
itok = input('');
if isempty(itok) || ~ismember(itok, 1:n)
  itok = 1; 
end
addpath([GSROOT '/tokamaks'])


% addpaths
tokpath = [d(itok).folder '/' d(itok).name];
addpath(tokpath)
addpath([tokpath '/eq'])
addpath([tokpath '/tok'])
addpath([tokpath '/matlab'])

% clear vars
clear tokpath d itok n i 







