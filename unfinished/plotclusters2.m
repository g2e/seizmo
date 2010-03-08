function [fh,sh]=plotclusters(data,grp,varargin)
%PLOTCLUSTERS    Plot clustered SEIZMO records
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar.  3, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  3, 2010 at 03:50 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:plotclusters:badNumInputs',...
        'Bad number of inputs!');
end

% deal with plotting options
P=plotparameters(varargin{:});

% select/open plot
if(isempty(P.FIGHANDLE) || P.FIGHANDLE<1)
    fh=figure;
else
    fh=figure(P.FIGHANDLE);
end

end
