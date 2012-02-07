function [varargout]=mt_g2c(g)
%MT_G2C    Extracts moment tensor components from the 3x3xN array format
%
%    Usage:    [M11,M22,M33,M12,M13,M23]=mt_g2c(momten)
%
%    Description:
%     [M11,M22,M33,M12,M13,M23]=MT_G2C(MOMTEN) extracts the components from
%     a 3x3xN moment tensor array as 6 Nx1 column vectors.  The output is
%     useful for expressions that deal with moment tensor components.
%
%    Notes:
%     - Throws an error if the tensor is not symmetric.
%
%    Examples:
%     % Convert from Harvard to Aki & Richards explicitly:
%     [mrr,mtt,mpp,mrt,mrp,mtp]=mt_g2c(momten);
%     [mxx,myy,mzz,mxy,mxz,myz]=deal(mtt,mpp,mrr,-mtp,mrt,-mrp)
%
%     % and check against the function that does this for you:
%     [mxx,myy,mzz,mxy,mxz,myz]=hrv2ar(mrr,mtt,mpp,mrt,mrp,mtp)
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_V2G, MT_G2V, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        June  1, 2011 - initial version
%        Feb.  7, 2012 - error if asymmetric, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 15:05 GMT

if(nargin==1)
    sz=size(g);
    if(~isreal(g) || ~isequal(sz(1:2),[3 3]))
        error('seizmo:mt_g2c:badInput',...
            'Input must be a real-valued 3x3xN array!');
    elseif(~isequal(g,permute(g,[2 1 3])))
        error('seizmo:mt_g2c:badInput',...
            'Cannot convert asymmetric tensors!');
    end
    g=permute(g,[3 1 2]);
    varargout={g(:,1) g(:,5) g(:,9) g(:,2) g(:,3) g(:,6)};
else
    error('seizmo:mt_g2c:badNumInput',...
        'Incorrect number of inputs!');
end

end