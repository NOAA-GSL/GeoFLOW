function [nl, Lmin, Lelem] = nicos(radius, res_or_ne, rnForm, p, q)
%
% function [nl, Lmin, Lelem]  = nicos(radius, res_or_ne, rnForm, p, q)
%
% Computes the number of 'Lagrange' divisions, nl, required
% for a given number of elements on an icos grid.
%
%  Usage:
%      [nl, Lmin, Lelem] = nicos(6738, 1.0, 'res');
%
%  Inputs:
%
%  radius   : sphere radius
%  res_or_ne: resolution (same units as radius) or number of elems,
%             depending on whether rnForm='res' or 'nelems'
%  rnForm   : flat telling what res_or_ne represents: if 'res', then
%             res_or_ne refers to the resolution; if 'nelems', then
%             res_or_ne is the number of elements.
%  p        : polynomial expansion order
%  q        : 'clustering' exponent, generally 2, but 1 <= q <= 2

%  Outputs:
%
%  nl       : number 'Lagrangian' division required
%  Lmin     : minimum length relative to radius resolvable
%  Lelem    : average element length

  if nargin < 1
    radius    = 1.0;
    res_or_ne = 0.01;
    rnForm    = 'res';
    p         = 1.0;
    q         = 2.0;
  end

  if nargin < 2
    res_or_ne = 0.01;
    rnForm    = 'res';
    p         = 1.0;
    q         = 2.0;
  end

  if nargin < 3
    rnForm    = 'res';
    p         = 1.0;
    q         = 2.0;
  end

  if nargin < 4
    p         = 1.0;
    q         = 2.0;
  end

  if nargin < 5
    q         = 2.0;
  end

  if ~(strcmp(lower(rnForm),'res') == 1 || strcmp(lower(rnForm),'nelems') == 1)
    error('Invalid value for rnForm');
  end

  % Generally, clustering exponent, q=2 describes te behavior
  % of Chebychev-like polynomials near the exdges of the
  % reference domain. p=1 would describe evenly spaced
  % nodes. We allow for any factor between these 2 extremes:
  if q < 1 || q > 2
    error('Invalid expansion exponent');
  end

  % For icos grid,
  % (1) Np = 10 nl^2 + 20 nl + 12, is number of triangle points given # Lagrange points
  % (2) Nt = 2*(Np-2) is the # triangles on grid
  % (3) Ne = 6(Np-2), is the number of elements given by 3 X the number of triangles
  
  stag = rnForm(find(~isspace(rnForm)))

  Area  = 4*pi*radius^2;  % global surf area
  if strcmp(lower(stag),'res') == 1   %  resolution specified
    Lmin  = res_or_ne;      % min len or res
    Lelem = p^q * Lmin;     % element length
    Ne    = Area / Lelem^2; % # elements required to cover surf evenly
    Np    = Ne/6 + 2;       % # triangle points from (3)
  elseif strcmp(lower(stag),'nelems') == 1 % nelems specified
    Ne    = res_or_ne;      % num elems
    Ae    = Area / Ne;      % avg area per elem
    Lelem = sqrt(Ae);       % avg elem length
    Lmin  = Lelem/(p^q);    % min len or res 
    Np    = Ne/6 + 2;       % # triangle points from (3)
  end

  % Solve (1) for # Lagrange partitions:
  a = 10;
  b = 20;
  c = 12 - Np;

  nl = (-b + sqrt(b^2 -4*a*c) ) / (2*a);

