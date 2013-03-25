function [strike,dip,rake]=normslip2sdr(normal,slip)
%NORMSLIP2SDR    Convert normal & slip vectors to strike, dip & rake angles
%
%    Usage:    sdr=normslip2sdr(normal,slip)
%              [strike,dip,rake]=normslip2sdr(normal,slip)
%
%    Description:
%     SDR=NORMSLIP2SDR(NORMAL,SLIP) converts normal and slip vectors to
%     strike, dip and rake angles.  NORMAL & SLIP are Nx3 arrays formatted
%     as [North East Up] where N allows for multiple vector sets to be
%     processed simultaneously.  NORMAL is the normal vector to the fault
%     and must point to the hanging wall (thus NORMAL(:,3) is always >=0).
%     The slip vector SLIP indicates the movement of the hanging wall
%     relative to the foot wall.  The output SDR is a Nx3 array arranged as
%     [strike dip rake] where strike is clockwise from North, dip is
%     positive downward from the horizontal and rake is counter-clockwise
%     in the fault plane from the strike direction.  Note that the strike
%     must be such that when you look along the direction of the strike the
%     fault dips to your right.
%
%     [STRIKE,DIP,RAKE]=NORMSLIP2SDR(NORMAL,SLIP) returns strike, dip &
%     rake separately.
%
%    Notes:
%
%    Examples:
%     % Vertical normal with a Northward slip:
%     normslip2sdr([0 0 1],[1 0 0])
%
%    See also: AUXPLANE, STRIKEDIP2NORM, NORM2STRIKEDIP, SDR2SLIP,
%              SDR2NULL, SDR2TPB, TPB2SDR, NODALLINES

%     Version History:
%        Mar. 15, 2013 - initial version
%        Mar. 18, 2013 - slip vector must be normalized, normal vector must
%                        point upwards
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% checks
sz1=size(normal,1); sz2=size(slip,1);
if(size(normal,2)~=3 || ndims(normal)>2)
    error('seizmo:normslip2sdr:badInput',...
        'NORMAL must be a Nx3 array as [North East Up] !');
elseif(~isnumeric(normal) || ~isreal(normal))
    error('seizmo:normslip2sdr:badInput',...
        'NORMAL must be a real-valued Nx3 array!');
elseif(size(slip,2)~=3 || ndims(slip)>2)
    error('seizmo:normslip2sdr:badInput',...
        'SLIP must be a Nx3 array as [North East Up] !');
elseif(~isnumeric(slip) || ~isreal(slip))
    error('seizmo:normslip2sdr:badInput',...
        'SLIP must be a real-valued Nx3 array!');
elseif(sz1~=sz2 && (sz1~=1 || sz2~=1))
    error('seizmo:normslip2sdr:badInput',...
        'NORMAL & SLIP different sizes, must be 1x3 or Nx3 arrays!');
end

% expand as needed
if(sz1==1); normal=normal(ones(sz2,1),:); end
if(sz2==1); slip=slip(ones(sz1,1),:); end

% check normal is upwards
if(any(normal(:,3)<0))
    error('seizmo:normslip2sdr:badNormal',...
        'Some normal vectors point downwards!');
end

% normalize slip vector
smag=sqrt(sum(slip.^2,2));
slip=slip./[smag smag smag];

% get strike & dip from normal
[strike,dip]=norm2strikedip(normal);

% get rake angle from slip vector and horizontal in-plane vector
rake=real(acosd(cosd(strike).*slip(:,1)+sind(strike).*slip(:,2)));

% fix sense of orientation as the above only returns positive angles
j=find(slip(:,3)<0);
rake(j)=-rake(j);

% combine if only one output
if(nargout<=1); strike=[strike(:) dip(:) rake(:)]; end

end
