function [nn,ne,en,ee]=split_horz_correlations(xc)
%SPLIT_HORZ_CORRELATIONS    Split horizontal correlations into nn/ne/en/ee
%
%    Usage:    [nn,ne,en,ee]=split_horz_correlations(xc)
%
%    Description:
%     [NN,NE,EN,EE]=SPLIT_HORZ_CORRELATIONS(XC) splits horizontal
%     correlations into NN, NE, EN, & EE sets where N & E stand for North &
%     East.  All other correlations are dropped.
%
%    Notes:
%     - This function is deprecated!  Please use HORZ_CORRELATIONS_SETS.
%     - The ordering nor the size of NN, NE, EN, & EE are not guaranteed to
%       match unless the original data records are ordered in a particular
%       way.  Using ROTATE immediately before CORRELATE is one way to
%       assure this (see the example below).
%
%    Examples:
%     % Read, Rotate, Correlate, Separate, Rotate:
%     data=readseizmo('*');
%     data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
%     xc=correlate(data,'mcxc');
%     [axc,xc]=split_auto_correlations(xc,false);
%     [nn,ne,en,ee]=split_horz_correlations(xc);
%     [rr,rt,tr,tt]=rotate_correlations(nn,ne,en,ee);
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS,
%              SPLIT_AUTO_CORRELATIONS

%     Version History:
%        Oct. 21, 2010 - initial version
%        Jan. 28, 2013 - doc update
%        Sep. 23, 2013 - fixed example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 23, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% only need orientation info
[scmp,mcmpinc,mcmpaz]=getheader(xc,'cmp','user2','user3');
nn=xc(mcmpinc==90 & mcmpaz==0  & scmp(:,1)==90 & scmp(:,2)==0 );
ne=xc(mcmpinc==90 & mcmpaz==0  & scmp(:,1)==90 & scmp(:,2)==90);
en=xc(mcmpinc==90 & mcmpaz==90 & scmp(:,1)==90 & scmp(:,2)==0 );
ee=xc(mcmpinc==90 & mcmpaz==90 & scmp(:,1)==90 & scmp(:,2)==90);

end

