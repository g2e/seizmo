function [varargout]=mt_v26(varargin)
%MT_V26    Converts moment tensor from Nx6 array to 6 Nx1 vectors
%
%    Usage:    [M11,M22,M33,M12,M13,M23]=mt_v26(mt)
%
%    Description: [M11,M22,M33,M12,M13,M23]=MT_V26(MT) converts a
%     moment tensor given as a Nx6 array to a separated format of Nx1
%     vectors.  The output format is useful for expressions that deal with
%     moment tensor components.
%
%    Notes:
%
%    Examples:
%     Convert from Harvard to Aki & Richards explicitly:
%      [mrr,mtt,mpp,mrt,mrp,mtp]=mt_v26(mt);
%      [mxx,myy,mzz,mxy,mxz,myz]=deal(mtt,mpp,mrr,-mtp,mrt,-mrp)
%     and check against the function that does this for you:
%      [mxx,myy,mzz,mxy,mxz,myz]=hrv2ar(mrr,mtt,mpp,mrt,mrp,mtp)
%
%    See also: MT_62V, MT_V2G, MT_G2V

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 21, 2010 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 21, 2010 at 11:25 GMT

% todo:

if(nargin==1)
    sz=size(varargin{1});
    if(~isreal(varargin{1}) || sz(2)~=6)
        error('seizmo:mt_v26:badInput',...
            'Input must be a real-valued Nx6 array!');
    end
    varargout=num2cell(varargin{1},1);
else
    error('seizmo:mt_v26:badNumInput',...
        'Incorrect number of inputs!');
end

end
