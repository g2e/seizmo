function [hd]=mo2hd(mo)
%MO2HD    Returns GlobalCMT half-duration estimate based on scalar moment
%
%    Usage:    hd=mo2hd(mo)
%
%    Description:
%     HD=MO2HD(MO) calculates the half-duration HD from the scalar moments
%     given in MO using the empirical relationship employed by the
%     GlobalCMT project since 1993.  This is great for getting half-
%     duration estimates for other catalogs (i.e. USGS MT).
%
%    Notes:
%     - Half-Duration scales an order of magnitude
%       per 3 orders in moment such that:
%        1e24 dyne-cm == 1.05s
%        1e27 dyne-cm == 10.5s
%
%    Examples:
%     % Here we compare entire catalog of GlobalCMT half-durations to their
%     % empirical scale.  Note that a significant portion of the catalog
%     % prior to 1993 has published half-durations following another trend
%     % (1e24 dyne-cm == 1.7s, 1e27 dyne-cm == 17s) giving significantly
%     % larger half-durations:
%     ocmts=findcmts('st',[1960 1],'et',[1993 1]);
%     ncmts=findcmts('st',[1993 1],'et',[2012 1]);
%     figure;
%     h1=subplot(2,1,1);
%     plot(h1,ocmts.scalarmoment.*10.^ocmts.exponent,ocmts.srcfuncdur ...
%         -mo2hd(ocmts.scalarmoment.*10.^ocmts.exponent),'.')
%     xlabel(h1,'scalar moment')
%     ylabel(h1,'published - empirical halfduration')
%     title(h1,'Pre-1993')
%     set(h1,'xscale','log');
%     h2=subplot(2,1,2);
%     plot(h2,ncmts.scalarmoment.*10.^ncmts.exponent,ncmts.srcfuncdur ...
%         -mo2hd(ncmts.scalarmoment.*10.^ncmts.exponent),'.')
%     xlabel(h2,'scalar moment')
%     ylabel(h2,'published - empirical halfduration')
%     title(h2,'1993-Present')
%     set(h2,'xscale','log');
%     linkaxes([h1 h2],'xy');
%     grid(h1,'on');
%     grid(h2,'on');
%
%    See also: READ_USGS_MT, FINDCMT, FINDCMTS

%     Version History:
%        June 10, 2011 - initial version
%        Mar. 19, 2013 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 19, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check scalarmoments
if(~isnumeric(mo) || ~isreal(mo) || any(mo<0))
    error('seizmo:mo2hd:badInput',...
        'MO must be an array of positive real values!');
end

% get half duration
%hd=1.7*(mo/1e24).^(1/3); % old GlobalCMT trend (pre-1993)
hd=1.05*(mo/1e24).^(1/3); % new GlobalCMT trend (1993-present)

end
