function [strike,dip,rake]=tpb2sdr(t,p,varargin)
%TPB2SDR    Returns the strike, dip & rake given principal axes
%
%    Usage:    sdr=tpb2sdr(t,p,b)
%              [strike,dip,rake]=tpb2sdr(t,p,b)
%
%    Description:
%     SDR=TPB2SDR(T,P,B) converts principal axes (tension, compression, &
%     null) to strike, dip, & rake of a focal mechanism.  T, P, & B must be
%     Nx3 arrays formatted as [value plunge azimuth] where N allows for
%     multiple focal mechanisms to be converted simultaneously.  Plunge and
%     azimuth are in degrees where plunge is positive downward from the
%     horizontal and azimuth is positive clockwise from North.  The output
%     SDR is a Nx3 array arranged as [strike dip rake] where strike is
%     clockwise from North, dip is positive downward from the horizontal
%     and rake is counter-clockwise in the fault plane from the strike
%     direction.  Note that the strike must be such that when you look
%     along the direction of the strike the fault dips to your right.
%
%     [STRIKE,DIP,RAKE]=TPB2SDR(T,P,B) returns strike, dip & rake
%     separately.
%
%    Notes:
%     - The null axis is not required and is an optional input.
%
%    Examples:
%     % What kind of focal mechanism makes
%     % maximum Pwaves directly North and East?
%     tpb2sdr([1 0 0],[1 0 90])
%
%    See also: SDR2TPB, STRIKEDIP2NORM, NORM2STRIKEDIP, SDR2SLIP, SDR2NULL,
%              NORMSLIP2SDR, NEU2VPA, VPA2NEU, NODALLINES, MT2TPB, TPB2MT

%     Version History:
%        Mar. 18, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% checking
sz1=size(t,1); sz2=size(p,1);
if(~isnumeric(t) || ~isreal(t) || ~isnumeric(p) || ~isreal(p))
    error('seizmo:tpb2sdr:badInput',...
        'T/P/B must be real-valued arrays!');
elseif(size(t,2)~=3 || ndims(t)>2 || size(p,2)~=3 || ndims(p)>2)
    error('seizmo:tpb2sdr:badInput',...
        'T/P/B must be Nx3 arrays as [VALUE PLUNGE AZIMUTH]!');
elseif(sz1~=sz2 && (sz1~=1 || sz2~=1))
    error('seizmo:tpb2sdr:badInput',...
        'T/P/B different sizes, must be 1x3 or Nx3 arrays!');
end

% expand as needed
if(sz1==1); t=t(ones(sz2,1),:); end
if(sz2==1); p=p(ones(sz1,1),:); end

% force to unit vectors
t(:,1)=1;
p(:,1)=1;

% convert to vector format (north, east, up)
t=vpa2neu(t);
p=vpa2neu(p);

% only need to get normal & slip
normal=t+p;
slip=t-p;

% precision fix
% - Converting between t & p axes of vertical faults (dip=90deg) and the
%   fault normal can give a slightly downwards normal.  We "fix" this by
%   adding a small amount to the vertical component.
fix=normal(:,3)<0 & normal(:,3)>-eps;
normal(fix,3)=normal(fix,3)+eps;

% if normal is pointing downwards, flip the sign of both the normal and
% slip vectors to give the appropriate normal and slip vectors (those
% corresponding to the hanging wall)
fix=normal(:,3)<0;
normal(fix,:)=-1*normal(fix,:);
slip(fix,:)=-1*slip(fix,:);

% get strike, dip & rake from normal & slip
[strike,dip,rake]=normslip2sdr(normal,slip);

% combine if only one output
if(nargout<=1); strike=[strike(:) dip(:) rake(:)]; end

end
