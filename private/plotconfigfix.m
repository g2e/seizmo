function [conf]=plotconfigfix(conf)
%PLOTCONFIGFIX    Fixes SEIZMO plot configuration control
%
%    Usage:    conf=plotconfigfix(conf)
%
%    Description: PLOTCONFIGFIX(CONF) fixes plot configuration structure
%     returned by PLOTCONFIG.  Currently this involves stepping through
%     a heirarchy of configuration fields and applying their values to
%     unset fields that are under them.  The heirarchy is actually defined
%     by PLOTCONFIG, so see that function for the layout.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     PLOTCONFIGFIX allows for changes to a high level field to change a
%     number of fields for instance:
%      plot1(data,'fgcolor','w','bgcolor','k')
%     uses two options that will apply their values to nearly every element
%     in the figure.
%
%    See also: plotconfig, plot1, plot2, plotall, recordsection, plotdendro

%     Version History:
%        Mar.  7, 2008 - initial version
%        Nov. 13, 2008 - renamed from PCONFFIX to PLOTCONFFIX
%        Nov. 15, 2008 - renamed from PLOTCONFFIX to PLOTCONFIGFIX
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 21:15 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end;

% check input
if(~isstruct(conf) || ~isfield(conf,'GENERAL'))
    error('seizmo:plotconfigfix:badInput',...
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
