function [varargout]=mt_v2c(varargin)
%MT_V2C    Converts moment tensor from Nx6 array to 6 Nx1 component vectors
%
%    Usage:    [M11,M22,M33,M12,M13,M23]=mt_v2c(mt)
%
%    Description:
%     [M11,M22,M33,M12,M13,M23]=MT_V2C(MT) extracts the components from a
%     Nx6 moment tensor array as 6 Nx1 column vectors.  The output format
%     is useful for expressions that deal with moment tensor components.
%
%    Notes:
%
%    Examples:
%     % Convert from Harvard to Aki & Richards explicitly:
%     [mrr,mtt,mpp,mrt,mrp,mtp]=mt_v2c(mt);
%     [mxx,myy,mzz,mxy,mxz,myz]=deal(mtt,mpp,mrr,-mtp,mrt,-mrp)
%
%     % and check against the function that does this for you:
%     [mxx,myy,mzz,mxy,mxz,myz]=hrv2ar(mrr,mtt,mpp,mrt,mrp,mtp)
%
%    See also: MT_C2V, MT_C2G, MT_V2G, MT_G2V, MT_G2C, MT_S2C, MT_S2G,
%              MT_S2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%        June  1, 2011 - doc update, renamed from mt_v26 to mt_v2c
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2011 at 15:05 GMT

% todo:

if(nargin==1)
    sz=size(varargin{1});
    if(~isreal(varargin{1}) || sz(2)~=6)
        error('seizmo:mt_v2c:badInput',...
            'Input must be a real-valued Nx6 array!');
    end
    varargout=num2cell(varargin{1},1);
else
    error('seizmo:mt_v2c:badNumInput',...
        'Incorrect number of inputs!');
end

end
