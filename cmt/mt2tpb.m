function [t,p,b]=mt2tpb(varargin)
%MT2TPB    Returns the principal axes of moment tensors
%
%    Usage:    [t,p,b]=mt2tpb(mt)
%              [t,p,b]=mt2tpb(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%
%    Description:
%     [T,P,B]=MT2TPB(MT) returns the principal axes of the moment tensor(s)
%     in MT.  MT must be a scalar struct as output by FINDCMT/FINDCMTS, a
%     Nx6 array, or a 3x3xN array where N is the number of moment tensors
%     in MT.  T, P, & B (the tension, pressure & null axes) are Nx3 arrays
%     of [eigenvalue plunge azimuth] where plunge & azimuth are in degrees.
%     Plunge is positive downward from the horizontal and azimuth is
%     positive clockwise from North.  Moment tensors in MT are expected to
%     be in Harvard convention (Up, South, East).
%
%     [T,P,B]=MT2TPB(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) allows specifying the moment
%     tensor components individually (Harvard system only).
%
%    Notes:
%     - Checking against the entire GlobalCMT catalog finds the following:
%       - eigenvalues match to within +/-.005, some to +/-.015, rare above
%       - plunges match to within +/-1, some to +/-2, rare go to +/-25
%       - azimuth match to within +/-1, some to +/-10, rare go to +/-225
%       The large azimuth deviations appear to be from roundoff
%       inaccuracies in plunge near 0 & 90.  This does not affect the
%       recomputed moment tensor much so all is well.
%
%    Examples:
%     % Get the principal axes of a CLVD moment tensor:
%     [t,p,b]=mt2tpb(diag([2 -1 -1]))
%
%    See also: TPB2MT, MT_DIAG, MT_UNDIAG, MT_DECOMP

%     Version History:
%        June 11, 2011 - initial version
%        June 13, 2011 - works with forced positive plunge
%        Mar. 19, 2013 - minor doc update
%        Nar, 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check tensor format (force 3x3xN)
error(mt_check(varargin{:}));
mt=mt_change('g',varargin{:});
n=size(mt,3);

% diagonalize moment tensor first
% (diagonal is eigenvalues, rotation vectors are the eigenvectors)
[mt,vec]=mt_diag(mt);

% put eigenvalues into Nx3 matrix [p b t]
t=[squeeze(mt(1,1,:)) squeeze(mt(2,2,:)) squeeze(mt(3,3,:))];

% flatten eigenvector matrix (3D to 2D)
vec=vec(:,:);

% convert eigenvectors in Harvard format to plunge and azimuth
[b,p]=cart2sph(vec(2,:),-vec(3,:),vec(1,:));
b=reshape(b,[3 n])';
p=reshape(p,[3 n])';

% delete mt, vec
clear mt vec;

% get into degrees
r2d=180/pi;
b=b*r2d;
p=p*r2d;

% require plunge is always positive
% - rotate azimuth 180 if plunge is negative
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

% now put v,p,a into the correct locations in t,p,b
[t,p,b]=deal([t(:,3) p(:,3) b(:,3)],...
    [t(:,1) p(:,1) b(:,1)],[t(:,2) p(:,2) b(:,2)]);

end
