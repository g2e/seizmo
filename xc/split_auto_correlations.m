function [axc,xc]=split_auto_correlations(xc,cmpflag)
%SPLIT_AUTO_CORRELATIONS    Split auto correlations from cross correlations
%
%    Usage:    [axc,xc]=split_auto_correlations(xc)
%              [axc,xc]=split_auto_correlations(xc,cmpflag)
%
%    Description:
%     [AXC,XC]=SPLIT_AUTO_CORRELATIONS(XC) splits auto correlations from
%     cross correlations.  Detection is based on name info (see CORRELATE
%     for the pertinent fields).  This is a simple convenience function.
%
%     [AXC,XC]=SPLIT_AUTO_CORRELATIONS(XC,CMPFLAG) ignores component codes
%     when deciding autocorrelations when CMPFLAG=FALSE.  This is useful
%     for detecting "auto" correlations between components of the same
%     station (ie. between horizontals).  The default is CMPFLAG=FALSE.
%     CMPFLAG=TRUE requires component codes to be exactly the same for
%     autocorrelation detection.
%
%    Notes:
%
%    Examples:
%     % Split off auto correlations so cross correlations can be reversed:
%     xc=correlate(data,'mcxc');
%     [axc,xc]=split_auto_correlations(xc);
%     rxc=reverse_correlations(xc);
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS, ISXC,
%              NAME_CORRELATIONS, NO_REDUNDANT_CORRELATIONS,
%              HORZ_CORRELATIONS_SETS, IS_FULL_MATRIX_OF_CORRELATIONS

%     Version History:
%        Oct. 21, 2012 - initial version
%        Jan. 28, 2013 - doc update
%        Sep. 20, 2013 - updated See also section, dropped time requirement
%        May  30, 2014 - cmpflag default switched to false as this is more
%                        generally useful
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2014 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default cmpflag
if(nargin<2 || isempty(cmpflag)); cmpflag=false; end
if(~isscalar(cmpflag) || ~islogical(cmpflag))
    error('seizmo:split_auto_correlations:badInput',...
        'CMPFLAG must be TRUE or FALSE!');
end

% autoxc = kname match
[kname,kt0,kt1,kt2,kt3]=getheader(xc,'kname','kt0','kt1','kt2','kt3');
if(cmpflag)
    auto=find(all(strcmp(kname,[kt0 kt1 kt2 kt3]),2));
else
    auto=find(all(strcmp(kname(:,1:3),[kt0 kt1 kt2]),2));
end
axc=xc(auto);
xc(auto)=[];

end

