function h = sphere2mercplot2d(svar, tindex, lon0, bwire, varargin)
%
% Plots 2d GeoFLOW spherical surface data to a plane via
% a Mercator projection.o
%
% Grid type must be GE_2DEMBEDDED.
%
%  Usage:
%    h  = sphere2mercplot2d('u1',10)
%
%  Input:
%    svar    : prefix for field file. Required
%    tindex  : time index for output. Required
%    lon0    : reference longitude in degrees. Default is 0 degrees.
%    bwire   : if > 0, print wire frame only; if 0 
%              print color patches. Default is 0
%    varargin: to pass to quadmesh: e.g. to plot
%              wire mesh only and set to single color, 
%              set bwire=1, and call:
%
%              sphere2mercplot2d('u1', 1, 90, 1, 'edgecolor','b')
%                 or
%              sphere2mercplot2d('u1',1, 90, 0,'colorbarlims',[-2 2])
%
%  Output:
%    h       : plot handle

%

if nargin < 2
  error('must specify svar, tindex');
end 
if nargin < 3
  lon0  = 0;
  bwire = 0;
end 
if nargin < 4
  bwire = 0;
end 

lon0 = lon0 * pi/180.0;

vartmp = varargin;
vartmp
bcolorbarlims = 0;
n = length(varargin)
j = 1;
while j <= n
  if strcmpi(vartmp{j},'colorbarlims') == 1
    if n < j+1
      error(sprintf('colorbarlims variable requires array [a b] of limits'));
    end
    colorbarlims = vartmp{j+1}
    bcolorbarlims = 1
    varargin = { vartmp{1:j-1}, vartmp{j+2:end} }
    j = j + 1
  end
  j = j + 1
end

scoord = {'xgrid','ygrid' 'zgrid'};

[umin, umax] = gminmax_gio(svar, tindex, 8, 'ieee-le');

d = dir('xgrid.*');
ntasks = length(d);
if ntasks<= 0 
  error('Grid data missing or incomplete');
end


nverts = 4;

%hwait = waitbar(0, 'Please wait...');

for itask = 0:ntasks-1

  % Read node coords:
  for j=1:3
    fname = sprintf('%s.%05d.out', scoord{j}, itask)
    [x{j} dim nelems porder gtype icycle time] = rgeoflow(fname, 8, 'ieee-le');
    if ( dim ~= 2 )
      error('Grid dimension must be 2');
    end
    if ( gtype ~= 2 )
      error('Grid type must be GE_2DEMBEDDED');
    end
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
    xx  = x{1}(icurr:icurr+lelem-1);
    yy  = x{2}(icurr:icurr+lelem-1);
    zz  = x{3}(icurr:icurr+lelem-1);
    uu  = u   (icurr:icurr+lelem-1);
    r   = sqrt(xx.^2 + yy.^2 + zz.^2);
    lat = asin(zz./r);
    lon = atan2(yy, xx));

    pdorder = double(porder);

    % Do mercator projection here:
    xxx = lon - lon0; 
    yyy = log2( tan(lat) + sec(lat) );


    % Find indices defining quads between
    % node points:
    imat = zeros(int64(prod(pdorder)/nverts), nverts);
    n = 1;
    for k = 1:pdorder(2) % cycle over each quad
      for j = 1:pdorder(1)
        bb = [ j+(k-1)*NN(1) j+1+(k-1)*NN(1) j+1+k*NN(1) j+k*NN(1) ];   
        imat(n,:) = bb;
        n = n + 1;
      end
    end
    cu = uu(imat);
    cm = mean(cu,2);
    cf = (cm - umin) / (umax - umin + eps);
    cv = (uu - umin) / (umax - umin + eps);
    if bwire == 0 
      h = quadmesh(imat,xxx,yyy,uu,'FaceColor','interp');
      colorbar('vertical');
    else
      h = quadmesh(imat,xxx,yyy,varargin{:});
    end
    hold on;

%   set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
    if bwire == 0
      title(sprintf('%s t=%f', svar, time));
    end
    icurr = icurr + lelem ; 

  end % end, elem loop
  
% waitbar(itask/ntasks,hwait);

end % end, task loop
%close(hwait);

if bcolorbarlims > 0 
  caxis(colorbarlims);
end

umin
umax
view(30,30)
axis equal
axis tight
%camlight
%lighting phong
%lighting gouraud
