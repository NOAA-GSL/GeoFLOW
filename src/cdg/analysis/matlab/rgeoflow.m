function [data dim nelems porder gtype time ivers] = rgeoflow(filein, isz, sformat)
%
% Reads binary GeoFLOW data, and stores in local variable data.
%
%  Usage:
%    data = rgeoflow(filename, 0, 8, 'ieee-be');
%
%  Input:
%    filein  : input file to read. Required.
%    isz     : data size (in bytes: either 4 or 8, e.g.). Default is isz=8.
%    sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%              endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8. Default
%              is 'ieee-le'.
%
%  Output:
%    data    : field data 
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
if nargin == 2
  sformat = 'ieee-le';
  swarn = sprintf('sformat=%s', sformat);
  warning(swarn);
end
if nargout < 1
  error('Must provide at least the data output argument');
end
if nargout > 7
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

% Read header:
[pdim pnelems pporder pgtype ptime pvers] = hgeoflow(filein, isz, sformat);
skip = sizeof(pvers) + sizeof(pdim) + sizeof(pnelems) + pdim*sizeof(pporder) + sizeof(pgtype) + sizeof(ptime);

lun =fopen(filein,'r',sformat);
if  lun == -1
  error(['File ' filein ' cannot be opened for reading']);
end
[fn permission thismachineformat] = fopen(lun); %machine format is for machine that reads, not that wrote
if ~strcmp(permission,'r')
   error('Invalid file')
end

if pdim ~= 2 && pdim ~= 3 
  error(['Invalid data dimension: ' dim]);
end

fseek(lun, skip, 'bof');
NN     = double(pporder + 1);

tsize = double(pnelems);
for j=1:pdim
  tsize = tsize * NN(j);
end
data = fread(lun, tsize, zsize);

fclose(lun);
if isempty(data)
  error(['File ' filein ' corrupted']);
end

 % Do output if required:
if nargout >= 2
  dim = pdim;
end
if nargout >= 3
  nelems = pnelems;
end
if nargout >= 4
  porder = pporder;
end
if nargout >= 5
  gtype = pgtype;
end
if nargout >= 6
  time = ptime;
end
if nargout == 7
  ivers = pvers;
end
