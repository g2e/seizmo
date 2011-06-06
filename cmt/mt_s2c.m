function [varargout]=mt_s2c(cmt)
%MT_S2C    Extracts moment tensor components from GlobalCMT struct
%
%    Usage:    [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=mt_s2c(cmt)
%
%    Description:
%     [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=MT_S2C(CMT) extracts the moment tensor
%     components from the GlobalCMT struct CMT.  This is useful for
%     operations that deal with moment tensor components.
%
%    Notes:
%
%    Examples:
%     % Convert from Harvard to Aki & Richards explicitly:
%     [mrr,mtt,mpp,mrt,mrp,mtp]=mt_s2c(cmt);
%     [mxx,myy,mzz,mxy,mxz,myz]=deal(mtt,mpp,mrr,-mtp,mrt,-mrp)
%
%     % and check against the function that does this for you:
%     [mxx,myy,mzz,mxy,mxz,myz]=hrv2ar(mrr,mtt,mpp,mrt,mrp,mtp)
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_V2G, MT_G2V, MT_G2C, MT_S2G,
%              MT_S2V

%     Version History:
%        June  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  1, 2011 at 11:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
if(~isstruct(cmt) || ~isscalar(cmt) || any(~isfield(cmt,...
        {'mrr' 'mtt' 'mpp' 'mrt' 'mrp' 'mtp'})))
    error('seizmo:mt_s2c:badInput',...
        'Input must be a scalar struct in the GlobalCMT format!');
end

% extract
varargout={cmt.mrr cmt.mtt cmt.mpp cmt.mrt cmt.mrp cmt.mtp};

end
