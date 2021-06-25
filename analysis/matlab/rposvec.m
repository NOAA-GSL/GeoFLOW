function [pos,vec] = rposvec(filein, isz, sformat)
%
% Reads binary GeoFLOW 'gbdy' data, and stores in local variable data.
%
%  Usage:
%    [pos, vec] = rposvec(filename, 8, 'ieee-be');
%
%  Input:
%    filein  : input file to read. Required.
%    isz     : data size (in bytes: either 4 or 8, e.g.). Default is isz=8.
%    sformat : data format of file: 'ieee-be' or 'ieee-le' for big-endian or little
%              endian if isz=4, or 'ieee-be.l64', 'ieee-le.l64' if isz=8. Default
%              is 'ieee-le'.
%
%  Output:
%    pos     : N x dim array of Cartesian positions; dim = 2 or 3 
%              depending on number of coords, and N is length of 
%              each component;
%    vec     : N x dim  array of Cartensian vector components.
%
if nargin < 1
  error('Input file name must be specified');
end
if nargin < 2
  isz = 8;
  sformat = 'ieee-le';
end
if nargin < 3
  sformat = 'ieee-le';
end

if nargout < 1
  error('Must provide pos and vec output arguments');
end
if nargout > 2
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
lun =fopen(filein,'r',sformat);
if  lun == -1
  error(['File ' filein ' cannot be opened for reading']);
end
[fn permission thismachineformat] = fopen(lun); %machine format is for machine that reads, not that wrote
if ~strcmp(permission,'r')
   error('Invalid file')
end

% Read header:
ncoords = fread(lun, 1      , 'uint32'); % version number
nvec    = fread(lun, 1      , 'uint32'); % problem dimension

ncoords = uint32(ncoords);
nvec    = uint32(nvec);


if ncoords ~= 2 && ncoords ~= 3 
  error(['Invalid data dimension: ' dim]);
end

pos = zeros(nvec, ncoords);
vec = zeros(nvec, ncoords);
for j = 1:ncoords % pos'n comps
  [pos(:,j), count] = fread(lun, nvec, zsize);
count
end
for j = 1:ncoords % vector comps
  [vec(:,j), count] = fread(lun, nvec, zsize);
count
end

fclose(lun);

end
