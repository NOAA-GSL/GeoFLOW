function [dim nelems porder gtype icycle time multivar ivers keys skip] = hgeoflow(filein, isz, sformat, quiet)
%
% Reads header from binary GeoFLOW data file
%
%  Usage:
%    [dim nelems porder gtype icycle time multivar ivers keys] = hgeoflow(filename, 0, 'ieee-be');
%
%  Input:
%    filein  : input file to read. Required.
%    isz     : data size (in bytes: either 4 or 8, e.g.). Default is isz=8.
%    sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%              endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8. Default
%              is 'ieee-le'.
%    quiet : if > 0, don't print header info. Default is 0.
%
%  Output:
%    dim     : data dimension (2, 3)
%    nelems  : number elements
%    porder  : array of size dim with the polynomial orders
%    gtype   : grid type (of GeoFLOW type GElemType)
%    icycle  : time cycle stamp
%    time    : time stamp
%    multivar: multiple fields in file?
%    ivers   : version number
%    keys    : element ids/keys
%    skip    : total header size in bytes
%
if nargin < 1
  error('Input file name must be specified');
end
if nargin < 2
  isz = 8;
  sformat = 'ieee-le';
  quiet = 0;
end
if nargin < 3
  sformat = 'ieee-le';
  quiet = 0;
end
if nargin < 4
  quiet = 0;
end

swarn = sprintf('using isz=%d; sformat=%s', isz, sformat);
warning(swarn);

if nargout > 10
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
pvers   = fread(lun, 1      , 'uint32'); % version number
pdim    = fread(lun, 1      , 'uint32'); % problem dimension
pnelems = fread(lun, 1      , 'uint64'); % # elems
pporder = fread(lun, pdim   , 'uint32'); % expansion order in each dir
pgtype  = fread(lun, 1      , 'uint32'); % grid type
pcycle  = fread(lun, 1      , 'uint64'); % time cycle 
ptime   = fread(lun, 1      ,  zsize  ); % time stamp
pmvar   = fread(lun, 1      , 'uint32'); % multiple vars in file
pkeys   = fread(lun, pnelems, 'uint64'); % element ids/keys

nkeys = pnelems;

% Ensure header types have correct size:
pvers   = uint32(pvers);
pdim    = uint32(pdim);
pnelems = uint64(pnelems);
pporder = uint32(pporder);
pgtype  = uint32(pgtype);
pcycle  = uint64(pcycle);
if strcmp(ssize,'real*4' )
  ptime = single(ptime);
elseif strcmp(ssize,'real*8')
  ptime = double(ptime);
end
pmultivar = uint32(pmvar);
pkeys     = uint64(pkeys);


pskip = sizeof(pvers) + sizeof(pdim)     + sizeof(pnelems) + pdim*sizeof(pporder) ...
                      + sizeof(pgtype)   + sizeof(pcycle)  + sizeof(ptime) ...
                      + sizeof(pmultivar) ...
                      + double(pnelems)*sizeof(pkeys)

pformat = '  %s=';
for j=1:pdim
  pformat = strcat(pformat,' %d');
end

if quiet == 0
  disp(sprintf('header for file: %s', filein));
  disp(sprintf('  %s=%d', 'vers'      , pvers));
  disp(sprintf('  %s=%d', 'dim'       , pdim));
  disp(sprintf('  %s=%d', 'nelems'    , pnelems));
  disp(sprintf(pformat  , 'pporder'   , pporder));
  disp(sprintf('  %s=%d', 'grid_type' , pgtype));
  disp(sprintf('  %s=%d', 'time_cycle', pcycle));
  disp(sprintf('  %s=%f', 'time_stamp', ptime));
  disp(sprintf('  %s=%f', 'multivar'  , pmultivar));
end

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
  icycle = pcycle;
end
if nargout >= 6
  time = ptime;
end
if nargout >= 7
  multivar = pmultivar;
end
if nargout >= 8
  ivers = pvers;
end
if nargout >= 9
  keys = pkeys;
end
if nargout == 10
  skip = pskip;
end


end
