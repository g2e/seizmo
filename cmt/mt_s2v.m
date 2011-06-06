function [mt]=mt_s2v(cmt)
%MT_S2V    Converts a GlobalCMT struct to a Nx6 moment tensor array
%
%    Usage:    mt=mt_s2v(cmt)
%
%    Description:
%     MT=MT_S2V(CMT) returns an Nx6 array containing the moment tensors in
%     the GlobalCMT struct CMT.  See functions READNDK, FINDCMT, & FINDCMTS
%     for more info on the GlobalCMT struct format.
%
%    Notes:
%     - GlobalCMT moment tensors are in Harvard form
%
%    Examples:
%     % Extract the moment tensors from the GloabCMT catalog:
%     mt=mt_s2v(findcmts);
%
%    See also: MT_C2V, MT_C2G, MT_V2C, MT_V2G, MT_G2V, MT_G2C, MT_S2C,
%              MT_S2G

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
    error('seizmo:mt_s2v:badInput',...
        'Input must be a scalar struct in the GlobalCMT format!');
end

% convert
mt=[cmt.mrr cmt.mtt cmt.mpp cmt.mrt cmt.mrp cmt.mtp]; % Nx6

end
