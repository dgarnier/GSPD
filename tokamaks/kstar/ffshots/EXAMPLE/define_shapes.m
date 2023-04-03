function shapes = define_shapes(opts)
%
%  inputs: opts struct
%          if opts.plotit==true, makes a gui plot of the target shapes
% 
%  outputs: the shape struct that has the target shapes vs time 
%
%     (rb.Data, rb.Time) - R position of boundary target shape, 
%                          timebase for R position of boundary target shape
%
%     (zb.Data, zb.Time) - Z position of boundary target shape
%     (rx.Data, rx.Time) - R position of target x-point
%     (zx.Data, zx.Time) - Z position of target x-point
%     (rtouch.Data, rtouch.Time) - R position of target touch point
%     (ztouch.Data, ztouch.Time) - Z position of target touch point
%     (rbdef.Data, rbdef.Time)   - R of boundary defining point
%     (zbdef.Data, zbdef.Time)   - Z of boundary defining point
%
%     Note that each of these can use different timebases which are later
%     interpolated. Each quantity should be defined at each time, but will 
%     only enter the optimization depending on the optimization weights. 
%
%


% In this example I mostly replicate the shapes from shot 23436
% with only slight modifications

efits = load('efits23436.mat').efits;
tok = load('kstar_tok.mat').tok;

tshapes = [0.8 1.2 1.5 1.8 2 2.3 3 3 10.1]; % which times to grabs shapes from
t =       [0.8 1.2 1.5 1.8 2 2.3 3 8 10];       % target times

tefit = [efits(:).time];
[~,i] = min(abs(tefit(:) - tshapes));
efits = efits(i);


% copy shape from the efits
rb = {};
zb = {};
rx = [];
zx = [];
N = length(efits);

for i = 1:N  
  eq = efits(i);
  k = eq.rbbbs == 0;    % remove the zero padding
  rb{i} = eq.rbbbs(~k);
  zb{i} = eq.zbbbs(~k); 
  [rx(i), zx(i)] = isoflux_xpFinder(eq.rg, eq.zg, eq.psizr, 1.4, -0.9);
end



% this code snippet is just an example on how one could modify the shape
% we will increase the elongation of the last equilibrium shape slightly
i = N;
s = shape_params(rb{i}, zb{i});

s.elong = s.elong + 0.02;  % modify fields of s as desired
                           % (elong, aminor, triu, etc)  

[rb{i}, zb{i}] = shape_edit(rb{i}, zb{i}, s);


%% Map target shapes to boundary-control points
% At this point, each of the (rb{i}, zb{i}) define a shape. Now sort these
% and map them onto control points. 


% interpolate to finer boundary
warning('off', 'MATLAB:polyshape:repairedBySimplify');
for i = 1:length(rb)
  [rb{i}, zb{i}] = interparc(rb{i}, zb{i}, 200, 1, 0);  % interpolate
  P = polyshape(rb{i}, zb{i});
  [rc,zc] = centroid(P);
  [rb{i}, zb{i}] = sort_ccw(rb{i}, zb{i}, rc, zc);      % sort 
end
  

% define control segments
segopts.rc = mean(tok.rg) - 0.05; 
segopts.zc = 0;
segopts.a = 0.25;
segopts.b = 0.35;
segopts.plotit = 0;
segopts.seglength = 4;
segs = gensegs(40, segopts);


% find intersections of boundary with segments
rcp = [];
zcp = [];
ncp = length(rb);
for i = 1:ncp
  [rcp(i,:), zcp(i,:)] = seg_intersections(rb{i}, zb{i}, segs, 0);  
end

shapes.rb.Time = t;
shapes.rb.Data = rcp;

shapes.zb.Time = t;
shapes.zb.Data = zcp;


%% define x-points, touch points
rmin = min(tok.limdata(2,:));

shapes.rx.Time = t;
shapes.rx.Data = rx;
% shapes.rx.Data = ones(size(t)) * 1.41;

shapes.zx.Time = t;
shapes.zx.Data = zx;
% shapes.zx.Data = ones(size(t)) * -0.88;

shapes.rtouch.Time = t([1 end]);
shapes.rtouch.Data = [rmin rmin]';

shapes.ztouch.Time = t([1 end]);
shapes.ztouch.Data = [0 0]';

% shapes.rbdef.Time = t;
% shapes.rbdef.Data = rx(:);
% shapes.rbdef.Data(t<=2) = rmin;
% 
% shapes.zbdef.Time = t;
% shapes.zbdef.Data = zx(:);
% shapes.zbdef.Data(t<=2) = 0;

shapes.rbdef.Time = t;
% shapes.rbdef.Data = min(shapes.rb.Data');
shapes.rbdef.Data = ones(size(t)) * 1.265;
shapes.zbdef.Time = t;
shapes.zbdef.Data = zeros(size(t));


shapes = check_structts_dims(shapes); 

if opts.plotit
  summary_shape_plot(shapes, tok);

  plot_structts(shapes, fields(shapes), 3);
  sgtitle('Shape targets'); drawnow
end


















