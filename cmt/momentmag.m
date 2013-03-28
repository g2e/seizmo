function [mw]=momentmag(varargin)
%MOMENTMAG    Returns the moment magnitude for a moment tensor
%
%    Usage:    Mw=momentmag(mt)
%              Mw=momentmag(m1,m2,m3,m4,m5,m6)
%
%    Description:
%     Mw=MOMENTMAG(MT) calculates moment magnitude Mw for the moment
%     tensor(s) given in MT.  MT must be a GlobalCMT struct
%     (formatted as returned by FINDCMTS), a Nx6 array or a 3x3xN array.
%     Note that for the numeric arrays the tensors should include the
%     exponent.  MT is expected to be in units of dyne*cm NOT N*m!
%
%     Mw=MOMENTMAG(M1,M2,M3,M4,M5,M6) allows inputting the moment tensor
%     components individually (Harvard or Aki & Richards system).
%
%    Notes:
%
%    Examples:
%     % Show the global seismicity rate as a function of magnitude,
%     % evaluate the completeness of the GlobalCMT catalog using the
%     % Gutenburg-Richter Law and show the associated b-value:
%     cmts=findcmts;
%     t=[cmts.year cmts.month cmts.day cmts.hour ...
%         cmts.minute cmts.seconds+cmts.centroidtime];
%     dt=timediff(t(1,:),t(end,:))/(86400*365.25); % in years
%     Mw=(4:.1:10)';
%     count=histc(momentmag(cmts),3.95:.1:10.05);
%     count(end)=[]; % drop Mw=10.05 value
%     fh=figure;
%     ax=axes;
%     pts=plot(ax,Mw,count/dt,'ko');
%     hold(ax,'on');
%     grid(ax,'on');
%     set(ax,'yscale','log');
%     xlabel('Magnitude (Mw)');
%     ylabel('#/yr');
%     title('GlobalCMT catalog Magnitude Completeness');
%     in54to8=find(Mw>=5.4 & Mw<=8);
%     p=polyfit(Mw(in54to8),log10(count(in54to8)/dt),1);
%     lh=plot(ax,Mw,10.^(Mw*p(1)+p(2)),'r','linewidth',2);
%     th=text(7,10,['b-value=' num2str(-p(1))],'parent',ax);
%     hold(ax,'off');
%
%    See also: FINDCMTS, FINDCMT, SCALARMOMENT

%     Version History:
%        Mar. 11, 2011 - initial version
%        Mar. 19, 2013 - doc update
%        Mar. 25, 2013 - update for mt_check/mt_change, minor fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check inputs
error(mt_check(varargin{:}));
if(isstruct(varargin{1})) % global cmt struct
    % use scalarmoment field as that is more accurate than using
    % the truncated moment tensor values which are always biased
    % to lower magnitudes
    mw=(2/3).*(log10(varargin{1}.scalarmoment ...
        .*10.^varargin{1}.exponent)-16.1);
else
    mw=(2/3).*(log10(scalarmoment(varargin{:}))-16.1);
end

end
