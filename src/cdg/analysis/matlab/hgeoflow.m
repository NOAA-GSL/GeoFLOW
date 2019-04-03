function [dim nelems porder gtype time ivers] = hgeoflow(filein, isz, sformat)
%
% Reads header from binary GeoFLOW data file
%
%  Usage:
%    [dim nelems porder gtype time ivers] = hghost(filename, 0, 'ieee-be');
%
%  Input:
%    filein  : input file to read. Required.
%    isz     : data size (in bytes: either 4 or 8, e.g.). Default is isz=8.
%    sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%              endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8. Default
%              is 'ieee-le'.
%
%  Output:
%    dim     : data dimension (2, 3)
%    nelems  : number elements
%    porder  : array of size dim with the polynomial orders
%    gtype   : grid type (of GeoFLOW type GElemType)
%    time    : time stamp
%    ivers   : version number
%
if nargin < 1
  error('Input file name must be specified');
end
if nargin == 1
  isz = 8;
  sformat = 'ieee-le';
  swarn = sprintf('using isz=%d; sformat=%s', isz, sformat);
  warning(swarn);
end
if nargout > 6
  error('Too many output arguments provided');
end

ssize = sprintf('real*%d',isz);
if strcmp(ssize,'real*4' )
  zsize = 'single';
elseif strcmp(ssize,'real*8')
  zsize = 'double';
else
  error('Type must be "real*4" or "real*8"');
end

%sformat
lun =fopen(filein,'r',sformat);
if  lun == -1
  error(['File ' filein ' cannot be opened for reading']);
end
[fn permission thismachineformat] = fopen(lun); %machine format is for machine that reads, not that wrote
if ~strcmp(permission,'r')
   error('Invalid file')
end

% Read header:
pvers   = fread(lun, 1   , 'uint32'); % version number
pdim    = fread(lun, 1   , 'uint32'); % problem dimension
pnelems = fread(lun, 1   , 'uint64'); % # elems
pporder = fread(lun, pdim, 'uint32'); % expansion order in each dir
pgtype  = fread(lun, 1   , 'uint32'); % grid type
ptime   = fread(lun, 1   ,  zsize  ); % time stamp

% Ensure header types have correct size:
pdim    = uint32(pdim);
pnelems = uint64(pnelems);
pporder = uint32(pporder);
pgtype  = uint32(pgtype);
if strcmp(ssize,'real*4' )
  ptime = single(ptime);
elseif strcmp(ssize,'real*8')
  ptime = double(ptime);
end

pformat = '  %s=';
for j=1:pdim
  pformat = strcat(pformat,' %d');
end

disp(sprintf('header for file: %s', filein));
disp(sprintf('  %s=%d', 'dim'       , pdim));
disp(sprintf('  %s=%d', 'nelems'    , pnelems));
disp(sprintf(pformat  , 'pporder'   , pporder));
disp(sprintf('  %s=%d', 'grid_type' , pgtype));
disp(sprintf('  %s=%d', 'time_stamp', ptime));

fclose(lun);

 % Do output if required:
if nargout >= 1
  dim = pdim;
end
if nargout >= 2
  nelems = pnelems;
end
if nargout >= 3
  porder = pporder;
end
if nargout >= 4
  gtype = pgtype;
end
if nargout >= 5
  ivers = pvers;
end
if nargout == 6
  time = ptime;
end
