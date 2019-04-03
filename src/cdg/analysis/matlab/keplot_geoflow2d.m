function h = keplot_geoflow2d(s1var, s2var, tindex)
%
% Does a mesh plot kinetic energy of 2D GeoFLOW data
%
%  Usage:
%    h = keplot_geoflow2d('u1','u2', 10)
%
%  Input:
%    s1var   : prefix for v1 component file. Required
%    s2var   : prefix for v2 component file. Required
%    tindex  : time index for output. Required
%
%  Output:
%    h       : plot handle

%

if nargin < 3 
  error('must specify s1var, s2var, and tindex');
end

scoord = {'xgrid','ygrid' 'zgrid'};

d = dir('xgrid.*');
ntasks = length(d);
if ntasks<= 0 
  error('Grid data missing or incomplete');
end

for itask = 0:ntasks-1

  % Read node coords:
  for j=1:2
    fname = sprintf('%s.%05d.out', scoord{j}, itask)
    [x{j} dim nelems porder gtype time] = rgeoflow(fname, 8, 'ieee-le');
  end


  fname = sprintf('%s.%06d.%05d.out', s1var, tindex, itask);
  [u1 dim nelems porder gtype time] = rgeoflow(fname, 8, 'ieee-le');
  fname = sprintf('%s.%06d.%05d.out', s2var, tindex, itask);
  [u2 dim nelems porder gtype time] = rgeoflow(fname, 8, 'ieee-le');
  if ( itask == 0 )
    figure;
  end

  NN = double(porder + 1);
  lelem = prod(NN);  % data length per element

  % Cycle over elems, and plot 'patches':
  icurr = 1;
  for n = 1:nelems
    xx = x{1}(icurr:icurr+lelem-1);
    yy = x{2}(icurr:icurr+lelem-1);
    uu1 = u1 (icurr:icurr+lelem-1);
    uu2 = u2 (icurr:icurr+lelem-1);
    xx = reshape(xx , NN(1), NN(2));
    yy = reshape(yy , NN(1), NN(2));
    uu1= reshape(uu1, NN(1), NN(2));
    uu2= reshape(uu2, NN(1), NN(2));
    surf( xx, yy, uu1.^2+uu2.^2 )
    title(sprintf('%s t=%f', 'KE', time));
    hold on
    icurr = icurr + lelem ; 
    
  end % end, elem loop
  

end % end, task loop

