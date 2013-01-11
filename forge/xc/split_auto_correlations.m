function [axc,xc]=split_auto_correlations(xc,cmpflag)
%SPLIT_AUTO_CORRELATIONS    Split auto correlations from cross correlations
%
%    Usage:    [axc,xc]=split_auto_correlations(xc)
%              [axc,xc]=split_auto_correlations(xc,cmpflag)
%
%    Description:
%     [AXC,XC]=SPLIT_AUTO_CORRELATIONS(XC) splits auto correlations from
%     cross correlations.  Detection is based on name info and start/end
%     times (see CORRELATE for the pertinent fields).  This is a simple
%     convenience function.
%
%     [AXC,XC]=SPLIT_AUTO_CORRELATIONS(XC,CMPFLAG) ignores component codes
%     when deciding autocorrelations when CMPFLAG=FALSE.  This is useful
%     for detecting "auto" correlations between components of the same
%     station.  The default is CMPFLAG=TRUE and requires component codes to
%     be exactly the same for autocorrelations.
%
%    Notes:
%
%    Examples:
%     % Split off auto correlations so cross correlations can be reversed:
%     xc=correlate(data,'mcxc');
%     [axc,xc]=split_auto_correlations(xc);
%     rxc=reverse_correlations(xc);
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS,
%              SPLIT_HORZ_CORRELATIONS

%     Version History:
%        Oct. 21, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 21, 2012 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default cmpflag
if(nargin<2 || isempty(cmpflag)); cmpflag=true; end
if(~isscalar(cmpflag) || ~islogical(cmpflag))
    error('seizmo:split_auto_correlations:badInput',...
        'CMPFLAG must be TRUE or FALSE!');
end

% autoxc = kname + abs times match
[kname,kt0,kt1,kt2,kt3,a,f,t0,t1]=getheader(xc,...
    'kname','kt0','kt1','kt2','kt3','a','f','t0','t1');
if(cmpflag)
    auto=find(all(strcmp(kname,[kt0 kt1 kt2 kt3]),2) & a==t0 & f==t1);
else
    auto=find(all(strcmp(kname(:,1:3),[kt0 kt1 kt2]),2) & a==t0 & f==t1);
end
axc=xc(auto);
xc(auto)=[];

end

