function [mt]=mt_s2g(cmt)
%MT_S2G    Converts a GlobalCMT struct to a 3x3xN moment tensor array
%
%    Usage:    momten=mt_s2g(cmt)
%
%    Description:
%     MOMTEN=MT_S2V(CMT) returns a 3x3xN array containing the moment
%     tensors in the GlobalCMT struct CMT.  See functions READNDK, FINDCMT,
%     & FINDCMTS for more info on the GlobalCMT struct format.
%
%    Notes:
%     - GlobalCMT moment tensors are in Harvard form
%
%    Examples:
%     % Extract the moment tensors from the GloabCMT catalog:
%     mt=mt_s2g(findcmts);
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_V2G, MT_G2V, MT_G2C, MT_S2C,
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
    error('seizmo:mt_s2g:badInput',...
        'Input must be a scalar struct in the GlobalCMT format!');
end

% convert
mt=permute(cat(3,[cmt.mrr cmt.mrt cmt.mrp],[cmt.mrt cmt.mtt cmt.mtp],...
    [cmt.mrp cmt.mtp cmt.mpp]),[2 3 1]); % 3x3xN

end
