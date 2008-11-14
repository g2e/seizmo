function [conf]=plotconffix(conf)
%PLOTCONFFIX    Fixes SAClab plot configuration control
%
%    Description: PLOTCONFFIX(CONF) fixes the plot configuration structure
%     returned by PLOTCONF.  Currently this just involves stepping through
%     a heirarchy of configuration fields and applying their values to
%     unset fields that are under them.  The heirarchy is actually defined
%     by PLOTCONF, so see that function for the layout.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    conf=plotconffix(conf)
%
%    Examples:
%     PLOTCONFFIX allows for changes to a high level field to change a
%     number of fields for instance:
%      plot1(data,'fgcolor','w','bgcolor','k')
%     uses two options that will apply their values to nearly every element
%     in the figure.
%
%    See also: plotconf, plot1, plot2, plotgrabbag, recsec, plotdendro

%     Version History:
%        Nov. 13, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2008 at 05:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check input
if(~isstruct(conf) || ~isfield(conf,'GENERAL'))
    error('SAClab:plotconffix:badInput',...
        'CONF must be a struct with a layout as from PLOTCONF!');
end

% heirarchy is stored in subfield GENERAL
fields=fieldnames(conf.GENERAL).';
for i=fields
    for j=conf.GENERAL.(i{:})
        if(isempty(conf.(j{:}))); conf.(j{:})=conf.(i{:}); end
    end
end

end
