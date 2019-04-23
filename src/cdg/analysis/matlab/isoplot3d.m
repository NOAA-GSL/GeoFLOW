function h = isoplot3(svar, tindex, isoval)
%
% Does a isosurface plot of 2D GeoFLOW data
%
%  Usage:
%    h  = isoplot3d('u1',10, 0.1)
%
%  Input:
%    s1var   : prefix for field file. Required
%    tindex  : time index for output. Required
%    isoval  : isosurface value
%
%  Output:
%    h       : plot handle

%

if nargin < 3
  error('must specify svar, tindex, isoval');
end 


ntasks = 2;
scoord = {'xgrid','ygrid' 'zgrid'};

d = dir('xgrid.*');
ntasks = length(d);
if ntasks<= 0 
  error('Grid data missing or incomplete');
end

umax = -1e35;
umin =  1e35;


for itask = 0:ntasks-1

  % Read node coords:
  for j=1:3
    fname = sprintf('%s.%05d.out', scoord{j}, itask)
    [x{j} dim nelems porder gtype icycle time] = rgeoflow(fname, 8, 'ieee-le');
  end
 

  fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  [u dim nelems porder gtype icycle time] = rgeoflow(fname, 8, 'ieee-le');
  if ( itask == 0 )
    figure;
  end

 
  NN = double(porder + 1);
  lelem = prod(NN(1:dim));  % data length per element

  % Cycle over elems, and plot 'patches':
  icurr = 1;
  for n = 1:nelems
    xx = x{1}(icurr:icurr+lelem-1);
    yy = x{2}(icurr:icurr+lelem-1);
    zz = x{3}(icurr:icurr+lelem-1);
    uu = u   (icurr:icurr+lelem-1);

    xx   = reshape(xx, NN(1), NN(2), NN(3));
    yy   = reshape(yy, NN(1), NN(2), NN(3));
    zz   = reshape(zz, NN(1), NN(2), NN(3));
    uu   = reshape(uu, NN(1), NN(2), NN(3));
    dx   = diff(xx,1,1); dxx = min(abs(dx(:)));
    dy   = diff(yy,1,2); dyy = min(abs(dy(:)));
    dz   = diff(zz,1,3); dzz = min(abs(dz(:)));

    xx   = reshape(xx, NN(1)* NN(2)* NN(3),1);
    yy   = reshape(yy, NN(1)* NN(2)* NN(3),1);
    zz   = reshape(zz, NN(1)* NN(2)* NN(3),1);
    uu   = reshape(uu, NN(1)* NN(2)* NN(3),1);
    pxyz = [xx yy zz];
    V    = uu;
    [pxyz, I] = unique(pxyz,'rows');
    V         = uu(I);

    umax = max(umax, min(abs(uu(:))));
    umin = min(umin, min(abs(uu(:))));

%sprintf('sort...')
    xm = unique(sort(xx(I))); 
    ym = unique(sort(yy(I))); 
    zm = unique(sort(zz(I)));
%sprintf('do meshgrid...')
    [xxi,yyi,zzi] = meshgrid(xm,ym,zm);
    pixyz = [xxi yyi zzi];
%sprintf('meshgrid done.')

if 0
    F   = scatteredInterpolant(pxyz, V);
%   F   = scatteredInterpolant(xx, yy, zz, uu);
    Vi   = F(xxi,yyi,zzi);
else
%sprintf('do interp...')
    F   = TriScatteredInterp(pxyz, V, 'natural');
%%  F   = TriScatteredInterp(xx, yy, zz, uu);
    Vi  = F(xxi,yyi,zzi);
%   Vi  = griddatan(pxyz, V, pixyz); 
%sprintf('interp done.')
end
%sprintf('do smoothing...')
%   Vi = smooth3(Vi, 'gaussian',5);
%sprintf('do isosurface...')
    p = patch( isosurface( xxi, yyi, zzi, Vi, isoval ) );
    isonormals(xxi, yyi, zzi, Vi, p);
%   p.FaceColor = 'blue';
%   p.EdgeColor = 'none';
    set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
    title(sprintf('%s t=%f: isoval=%f', svar, time, isoval));
    hold on
    icurr = icurr + lelem ; 

    
  end % end, elem loop
  

end % end, task loop

sprintf('data range: (%f, %f)', umin, umax)
if isoval > umax || isoval < umin 
  warning(sprintf('Desired isoval=%f ourside data range',isoval));
end

view(3)
axis equal
axis tight
camlight
lighting phong
%lighting gouraud
