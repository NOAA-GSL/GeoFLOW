function [c h levels umin umax] = boxcont2d(svar, tindex, levtype, xlev, blog, dtype, isz )
%
%   function h = boxcont2d(svar, tindex, blog, xlev, levtype, nozero, dtype, isz)
%
% Plots contour plots for 2D GeoFLOW Posix or collective data. 
%
%  Usage:
%    [c h] = boxcont2d('u1',10, 0, 8)
%
%  Input:
%    svar    : prefix for field file. Required
%    tindex  : time index for output. Required
%    levtype : if 'num', then xlev is the number of levels to use.
%              If levtype=='delta', then xlev is the interval
%              between contour levels. If levtype=='lev', then
%              xlev is the array of contour levels. The default
%              is levtype=='num', with 5 levels.
%    xlev    : number of levels, if levtype = 'num'; contour
%              increment if levtype = 'delta', actual contour levels
%              if levtype=='lev'. Default if levtype=='num'; else
%              must be specified.
%    blog    : take log of data? Default is 0
%    dtype   : data file type: 'POSIX', 'COLL'. Default is 'COLL'
%    isz    : floating point size (4 or 8). Default is 8.
%
%  Output:
%    c       : contour matrix
%    h       : plot handle
%
if nargin < 2
  error('must specify svar and tindex');
end 

if nargin < 3
  levtype = 'num';
  xlev    = 5;
  blog    = 0;
  dtype   = 'COLL';;
  isz     = 8;
end 
if nargin < 4
  if ~strcmpi(levtype, 'num')
    error(['xlev must be set when levtype=' levtype]);
  end
  xlev    = 5;
  blog    = 0;
  dtype   = 'COLL';;
  isz     = 8;
end 
if nargin < 5
  blog    = 0;
  dtype   = 'COLL';;
  isz     = 8;
end 
if nargin < 6
  dtype   = 'COLL';;
  isz     = 8;
end 
if nargin < 7
  isz     = 8;
end 
if nargin < 8
  isz = 8;
end 

if  ~strcmpi(levtype,'lev') && ~strcmpi(levtype,'delta') ...
 && ~strcmpi(levtype,'num')
    error(['levtype value ', levtype, ' invalid']);
end

dtype
if ~strcmp(dtype,'POSIX') & ~strcmp(dtype,'COLL')
  error(['Invalid dtype: ' dtype]);
end

if strcmpi(levtype, 'lev')
  levels = xlev;
end

lwidth = 2;
szfont = 16;

bcolorbarlims = 0;

scoord = {'xgrid','ygrid' 'zgrid'};

[umin, umax] = gminmax_gio(svar, tindex, dtype, isz, 'ieee-le');

ntasks = 1;
if strcmp(dtype,'POSIX')
  d = dir('xgrid.*');
  ntasks = length(d);
  if ntasks<= 0 
    error('Grid data missing or incomplete');
  end
end

nverts = 4;

figure();
for itask = 0:ntasks-1

  % Read node coords:
  for j=1:2
    if strcmp(dtype,'POSIX')
      fname = sprintf('%s.%06d.%05d.out', scoord{j}, 0, itask);
    elseif dtype == 'COLL'
      fname = sprintf('%s.%06d.out', scoord{j}, 0);
    end
    [x{j} dim nelems porder gtype icycle time mvar] = rgeoflow(fname, isz, 'ieee-le');
  end
  if ( dim ~= 2 )
    error('Grid must be 2D');
  end 

  if strcmp(dtype,'POSIX')
    fname = sprintf('%s.%06d.%05d.out', svar, tindex, itask);
  elseif dtype == 'COLL'
    fname = sprintf('%s.%06d.out', svar, tindex);
  end
  [u dim nelems porder gtype icycle time mvar] = rgeoflow(fname, isz, 'ieee-le');


 umax = max(u);
 umin = min(u);
 
  NN = double(porder + 1); 
  if blog == 1 
    u = log10(abs(u));
  end

if 0
   % print max(u) and its position:
   aumax = max(abs(umax),abs(umin));
   aumin = min(abs(umax),abs(umin));
   I = find(abs(u) > 0.95*aumax & abs(u) <= aumax);
   for j=1:numel(I)
     sprintf('(%f, %f); u=%f',x{1}(I(j)), x{2}(I(j)), u(I(j)) )
   end
end

  xorig = [x{1}(:) x{2}(:)];
  [xd, IA, IC]  = unique(xorig,'rows');
  clear x ;
  [x z] = deal(xd(:,1), xd(:,2));
  clear xd xorig;
  U = u(IA);
  clear u;

  dx = diff(x);
  I  = find(abs(dx) > 5.0e-4*abs(max(x)-min(x)) );
  ngridx = 1.25*( max(x) - min(x) ) / min(abs(dx(I)));
  ngridx = int32(ngridx);

  dx = diff(z);
  I  = find(abs(dx) > 5.0e-4*abs(max(x)-min(x)) );
  ngridy = 1.1*( max(z) - min(z) ) / min(abs(dx(I)));
  ngridy = int32(ngridy)
 
  [X,Y] = ndgrid(linspace(min(x),max(x),ngridx),linspace(min(z),max(z),ngridy));

  Z     = griddata(x,z,U(:),X,Y,'linear');

if 1
  % In case there's topography, we set values to nan 'inside' it:
  K       = boundary(x, z, 1.0);
  PGON    = polyshape(x(K), z(K),'Simplify', false);
%plot(PGON)
%whos PGON
%fieldnames(PGON,'-full')
  VERTS   = getfield(PGON,'Vertices');

% IDX     = isinterior(PGON, X(:), Y(:));
  IDX     = inpolygon(X,Y,VERTS(:,1),VERTS(:,2));
  Z(~IDX) = nan;

numel(IDX)
numel(Z)
end


  zmin = min(min(Z));
  zmax = max(max(Z));
  if strcmpi(levtype,'num')
    nlev = int32(xlev);
    dc = (zmax - zmin)/(double(nlev)+1);
    dcvec = [0:double(nlev)].*dc + zmin;
  elseif strcmpi(levtype,'delta')
    nlev = int32((zmax - zmin)/xlev);
    dc = (zmax - zmin)/(double(nlev)+1);
    dcvec = [0:double(nlev)].*dc + zmin;
  end

  % Remove zero-contour:
  if ~strcmpi(levtype,'lev')
    levels = dcvec;
  end

  [c h] = contourf(X, Y, Z, levels); 
  axis normal;

  hold on;


end % end, task loop

if 1
% Read normal vectors and plot:
[pos, vec] = rposvec('gbdy.out', isz)
dn = 1;
quiver(pos(1:dn:end,1), pos(1:dn:end,2), vec(1:dn:end,1), vec(1:dn:end,2), 0.1);
end

sminmax = sprintf('min=%0.5g; max=%0.5g', umin, umax);
stitle  = sprintf('%s, t=%0.5g: \n %s', svar, time, sminmax);
title(stitle);
xlabel('x (m)');
ylabel('z (m)');


end
