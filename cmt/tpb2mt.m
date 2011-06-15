function [mt]=tpb2mt(t,p,b)
%TPB2MT    Returns moment tensors given their principal axes
%
%    Usage:    mt=tpb2mt(t,p,b)
%
%    Description:
%     MT=TPB2MT(T,P,B) assembles the principal axes info to recreate the
%     associated moment tensor(s).  T, P, & B should be Nx3 arrays of
%     [eigval plunge azimuth] where plunge & azimuth are in degrees.  The
%     moment tensor is returned in Harvard convention.
%
%    Notes:
%
%    Examples:
%     % Compare T, P, & B info from GlobalCMT catalog & calculated:
%     cmt=findcmts;
%     t=[cmt.eigval1 cmt.plunge1 cmt.azimuth1];
%     p=[cmt.eigval3 cmt.plunge3 cmt.azimuth3];
%     b=[cmt.eigval2 cmt.plunge2 cmt.azimuth2];
%     cmt=mt_s2g(cmt);
%     max(max(max(cmt-tpb2mt(t,p,b))))
%     [t,p,b]=mt2tpb(cmt);
%     max(max(max(cmt-tpb2mt(t,p,b))))
%
%    See also: MT2TPB, MT_DIAG, MT_UNDIAG, MT_DECOMP

%     Version History:
%        June 12, 2011 - initial version
%        June 13, 2011 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2011 at 13:50 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% check principal axes
if(~isnumeric(t) || ~isreal(t) || size(t,2)~=3 || ndims(t)~=2)% ...
%        || any(t(:,2)<0 | t(:,2)>90))
    error('seizmo:tpb2mt:badInput','T must be a real-valued Nx3 array!');
elseif(~isnumeric(p) || ~isreal(p) || size(p,2)~=3 || ndims(p)~=2)% ...
%        || any(p(:,2)<0 | p(:,2)>90))
    error('seizmo:tpb2mt:badInput','P must be a real-valued Nx3 array!');
elseif(~isnumeric(b) || ~isreal(b) || size(b,2)~=3 || ndims(b)~=2)% ...
%        || any(b(:,2)<0 | b(:,2)>90))
    error('seizmo:tpb2mt:badInput','B must be a real-valued Nx3 array!');
elseif(~isequal(size(t),size(p),size(b)))
    error('seizmo:tpb2mt:badInput','T, P, & B must be equal sized!');
end
n=size(t,1);

% plunge, azimuth => eigenvectors
vec=nan(n,3,3); d2r=pi/180;
[vec(:,:,2),vec(:,:,3),vec(:,:,1)]=sph2cart(...
    -d2r*[p(:,3) b(:,3) t(:,3)],d2r*[p(:,2) b(:,2) t(:,2)],ones(n,3));
vec=permute(vec,[3 2 1]);

% make diagonal matrix from eigenvalues
mt=zeros(3,3,n);
mt(1,1,:)=reshape(p(:,1),[1 1 n]);
mt(2,2,:)=reshape(b(:,1),[1 1 n]);
mt(3,3,:)=reshape(t(:,1),[1 1 n]);

% eigenvalues & eigenvectors => mt
mt=mt_undiag(mt,vec);

end
