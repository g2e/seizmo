function h=m_patch(long,lat,C,varargin);
% M_PATCH Create patches on a map
%    M_PATCH(LONG,LAT,C) is a drop-in replacement for PATCH that uses 
%    longitude/latitude coordinates to draw a patch on the current map. 
%    See PATCH for more details about the way in which patch colours and 
%    properties should be specified.
%
%    Currently you cannot specify C to be other than a string or 1x3 RGB
%    vector.
%
%    See also M_LINE, M_LL2XY

% Rich Pawlowicz (rich@ocgy.ubc.ca) 3/Sep/98
% 
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

%  10/Mar/99 - changed order of calls ('c' not handled correctly in mu_coast otherwise)
% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
% 8/Feb/11 - draw at +/-360 too (gge)

% row vector to column vector
if(isvector(long)); long=long(:); end
if(isvector(lat)); lat=lat(:); end

% draw +/-360
long=[long long+360 long-360];
lat=[lat lat lat];
[m,n]=size(long);
h=mu_coast('vector',[reshape([long;long(1,:);NaN+ones(1,n)],(m+2)*n,1),...
    reshape([lat;lat(1,:);NaN+ones(1,n)],(m+2)*n,1)],'patch',C,...
    'tag','m_patch',varargin{:});

if nargout==0,
 clear h
end;
