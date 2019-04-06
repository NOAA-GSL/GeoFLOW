function h = keplot_geoglow2d(s1var, s2var, tibeg, tiend, iskip)
%
% Make frames of a movie of KE 
%
if nargin < 4
  error('must specify s1var, s2var, tibeg, and tiend');
end 
if nargin < 5
  iskip = 1;
end

for it = tibeg:iskip:tiend
  keplot_geoflow2d(s1var, s2var, it);
  fname = sprintf('%s.%06d.pdf', 'ke_per', it)
  print('-dpdf', fname);
end % end, time loop

