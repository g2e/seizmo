function [mt]=mt_undiag(mt,vec)
%MT_UNDIAG    De-diagonalizes moment tensor(s) into Harvard orientation
%
%    Usage:    mt=mt_undiag(val,vec)
%
%    Description:
%     MT=MT_UNDIAG(VAL,VEC) undoes the diagonalization (aka rotation) of
%     the moment tensor(s) in VAL using the eigenvectors in VEC.  This
%     rotates them from the principal axes basis back into the Harvard
%     convention (Up, South, East).  This is useful for plotting
%     decomposed moment tensors (deviatoric, best double couple,  clvd,
%     etc).  VAL & VEC should follow the convention returned by MT_DIAG,
%     which is that they must both be 3x3xN.  MT is 3x3xN.
%
%    Notes:
%
%    Examples:
%     % Decompose some cmts into their major & minor double-couples,
%     % undiagonalize and plot:
%     cmts=findcmt('n',10);
%     [major,minor,vec]=mt_decomp(cmts,'majmin');
%     major=mt_undiag(major,vec);
%     minor=mt_undiag(minor,vec);
%     figure;
%     plotmt(1:10,2,cmts)
%     hold on;
%     plotmt(1:10,1,major)
%     plotmt(1:10,0,minor)
%     axis equal tight off
%
%    See also: MT_DIAG, MT_DECOMP, MT2TPB, TPB2MT, PLOTMT

%     Version History:
%        June 11, 2011 - initial version
%        Mar. 20, 2013 - doc update, floating point asymmetry fix
%        Mar. 25, 2013 - fixed example for mt_change/mt_check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check eigenvalues & vectors
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt_undiag:badInput',...
        'VAL must be a 3x3xN real-valued matrix!');
elseif(~isnumeric(vec) || ~isreal(vec))
    error('seizmo:mt_undiag:badInput',...
        'VEC must be a 3x3xN real-valued matrix!');
elseif(~isequal(size(mt),size(vec)))
    error('seizmo:mt_undiag:badInput',...
        'VAL & VEC are not equal sized!');
elseif(size(mt,1)~=3 || size(mt,2)~=3 || ndims(mt)>3)
    error('seizmo:mt_undiag:badInput',...
        'VAL & VEC are not 3x3xN sized matrices!');
end

% convert eigenvectors & eigenvalues to harvard
for i=1:size(mt,3)
    mt(:,:,i)=vec(:,:,i)*mt(:,:,i)*vec(:,:,i)';
end

% floating point asymmetry fix
mt=(mt+permute(mt,[2 1 3]))/2;

end
