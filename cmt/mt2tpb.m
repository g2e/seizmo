function [t,p,b]=mt2tpb(mt)
%MT2TPB    Returns the principal axes of moment tensors
%
%    Usage: [t,p,b]=mt2tpb(mt)
%
%    Description:
%     [T,P,B]=MT2TPB(MT) returns the principal axes of the moment tensor(s)
%     in MT.  MT must be a Nx6 or 3x3xN array where N is the number of
%     moment tensors in MT.  T, P, & B are Nx3 arrays of [val plunge azi]
%     where plunge & azi are in degrees.  To recover the moment tensor see
%     TPB2MT.
%
%    Notes:
%     - Checking against the entire GlobalCMT catalog finds the following:
%       - eigenvalues match to within +/-.005, some to +/-.015, rare above
%       - plunges match to within +/-1, some to +/-2, rare go to +/-25
%       - azimuth match to within +/-1, some to +/-10, rare go to +/-225
%       The large azimuth deviations appear to be from roundoff
%       inaccuracies in plunge near 0 & 90.  This affects the recomputed
%       moment tensor little so all is well.
%
%    Examples:
%     % Get the principal axes of a clvd moment tensor:
%     [t,p,b]=mt2tpb(diag([2 -1 -1])
%
%    See also: TPB2MT, MT_DIAG, MT_UNDIAG, MT_DECOMP

%     Version History:
%        June 11, 2011 - initial version
%        June 13, 2011 - works with forced positive plunge
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2011 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check tensor format (force 3x3xN)
mtsz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt2tpb:badInput',...
        'MT must be a real-valued numeric array!');
elseif(isequal(mtsz(1:2),[3 3]) && any(numel(mtsz)==[2 3]))
    if(numel(mtsz)>2); n=mtsz(3); else n=1; end
elseif(mtsz(2)==6 && numel(mtsz)==2)
    mt=mt_v2g(mt); % convert from Nx6 to 3x3xN
    n=mtsz(1);
else
    error('seizmo:mt2tpb:badInput',...
        'MT must be a harvard moment tensor array as 3x3xN or Nx6!');
end

% diagonalize first
[mt,vec]=mt_diag(mt);

% put eigenvalues into Nx3 matrix [p b t]
t=[squeeze(mt(1,1,:)) squeeze(mt(2,2,:)) squeeze(mt(3,3,:))];

% flatten vec
vec=vec(:,:);

% convert eigenvectors to plunge and azimuth
[b,p]=cart2sph(vec(2,:),-vec(3,:),vec(1,:));
b=reshape(b,[3 n])';
p=reshape(p,[3 n])';

% delete mt, vec
clear mt vec;

% get into degrees
r2d=180/pi;
b=b*r2d;
p=p*r2d;

% plunge is always positive so rotate azimuth 180 if plunge is negative
b(p<0)=b(p<0)+180;
p=abs(p);

% force azimuth to be in 0-360
b(b>=360)=b(b>=360)-360;
b(b<0)=b(b<0)+360;

% recover eigenvectors (for debugging)
%[y,z,x]=sph2cart(-b/r2d,p/r2d,ones(size(b)));
%for i=1:n
%    [x(i,:); y(i,:); z(i,:)]
%end

% v,p,a => t,p,b
[t,p,b]=deal([t(:,3) p(:,3) b(:,3)],...
    [t(:,1) p(:,1) b(:,1)],[t(:,2) p(:,2) b(:,2)]);

end
