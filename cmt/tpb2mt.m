function [mt]=tpb2mt(t,p,b)
%TPB2MT    Returns moment tensors given their principal axes
%
%    Usage:    mt=tpb2mt(t,p,b)
%
%    Description:
%     MT=TPB2MT(T,P,B) assembles the principal axes info (tension,
%     pressure, and null axes) to recreate the associated moment tensor(s).
%     T, P, & B should be Nx3 arrays of [eigenvalue plunge azimuth] where
%     N allows for multiple moment tensors to be processed simultaneously.
%     Plunge & azimuth are in degrees where plunge is positive downward
%     from the horizontal and azimuth is positive clockwise from North.
%     The moment tensors are returned as a 3x3xN array in Harvard
%     convention (Up, South, East).
%
%    Notes:
%
%    Examples:
%     % Comparing differences in moment tensor components published in the
%     % GlobalCMT catalog & that calculated from their principal axes info:
%     cmts=findcmts;
%     t=[cmts.eigval1 cmts.plunge1 cmts.azimuth1];
%     p=[cmts.eigval3 cmts.plunge3 cmts.azimuth3];
%     b=[cmts.eigval2 cmts.plunge2 cmts.azimuth2];
%     max(mt_change('v',cmts)-mt_change('v',tpb2mt(t,p,b)))
%
%     % Now compare to the accuracy of the conversion back and forth:
%     [t,p,b]=mt2tpb(cmts);
%     max(mt_change('v',cmts)-mt_change('v',tpb2mt(t,p,b)))
%
%    See also: MT2TPB, MT_DIAG, MT_UNDIAG, MT_DECOMP

%     Version History:
%        June 12, 2011 - initial version
%        June 13, 2011 - added docs
%        Mar. 19, 2013 - doc update, minor fix for precision issue
%        Mar. 25, 2013 - fix examples for mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% check principal axes
if(~isnumeric(t) || ~isreal(t) || size(t,2)~=3 || ndims(t)~=2)
    error('seizmo:tpb2mt:badInput','T must be a real-valued Nx3 array!');
elseif(~isnumeric(p) || ~isreal(p) || size(p,2)~=3 || ndims(p)~=2)
    error('seizmo:tpb2mt:badInput','P must be a real-valued Nx3 array!');
elseif(~isnumeric(b) || ~isreal(b) || size(b,2)~=3 || ndims(b)~=2)
    error('seizmo:tpb2mt:badInput','B must be a real-valued Nx3 array!');
elseif(~isequal(size(t),size(p),size(b)))
    error('seizmo:tpb2mt:badInput','T, P, & B must be equal sized!');
end
n=size(t,1);

% plunge & azimuth => eigenvectors (in Harvard convention)
vec=nan(n,3,3); d2r=pi/180;
[vec(:,:,2),vec(:,:,3),vec(:,:,1)]=sph2cart(...
    -d2r*[p(:,3) b(:,3) t(:,3)],d2r*[p(:,2) b(:,2) t(:,2)],ones(n,3));
vec=permute(vec,[3 2 1]);

% make diagonal matrix from eigenvalues
mt=zeros(3,3,n);
mt(1,1,:)=reshape(p(:,1),[1 1 n]);
mt(2,2,:)=reshape(b(:,1),[1 1 n]);
mt(3,3,:)=reshape(t(:,1),[1 1 n]);

% eigenvalues & eigenvectors => mt (in Harvard convention)
mt=mt_undiag(mt,vec);

% floating point fix
% - tensor is not symmetric by a very small amount (<10^-14)
% - workaround this by averaging tensor w/ transpose
mt=(mt+permute(mt,[2 1 3]))/2;

end
