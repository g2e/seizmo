function [data]=autowindow(data,varargin)
%AUTOWINDOW    Automatic windowing of SEIZMO records
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: GETPEAKS, ENVELOPE

% todo:
% - what options are available:
%   - seedtime
%     - record max by default
%   - mindrop
%     - 0.25 by default
%   - maxwidth
%     - record limits by default
%   - maxbump
%     - 0.25 by default
%   - problem
%     - delete - delete records that fail mindrop
%     - width - return entire maxwidth for records that fail mindrop
%     - min - return window from minimums for records that fail mindrop

% - how do we go from seed time to peak time
%   - how do we decide if peak is just a minor peak
%     - if it cannot pass mindrop/maxbump in maxwidth
%       - 

end
