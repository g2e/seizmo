function [s,freq]=sea_spectrum(u10,freq,fetch,type)
%SEA_SPECTRUM    Returns a parametric 1D wave spectrum for sea waves
%
%    Usage:    [s,freq]=sea_spectrum(u10)
%              s=sea_spectrum(u10,freq)
%              s=sea_spectrum(u10,freq,fetch,'jonswap')
%
%    Description:
%     [S,FREQ]=SEA_SPECTRUM(U10) returns the Pierson-Moskowitz 1D
%     parametric sea spectrum for a given windspeed U10 in meters per
%     second (note that U10 is the windspeed measured 10 meters above the
%     surface).  S is a power spectral density in m^2/Hz and FREQ are the
%     associated frequencies in Hz.  FREQ is regularly sampled in base-10
%     logarithmic space over the frequency range 0.01-1Hz.  If U10 is a
%     vector then a spectrum is computed for each element in U10 and S(:,i)
%     gives the spectrum for U10(i).
%
%     S=SEA_SPECTRUM(U10,FREQ) calculates the Pierson-Moskowitz spectrum at
%     the given frequencies FREQ.
%
%     S=SEA_SPECTRUM(U10,FREQ,FETCH,'JONSWAP') calculates the JONSWAP 1D
%     parametric sea spectrum for a given windspeed & fetch FETCH in
%     kilometers.  FREQ may be given or set to [] to use the default.
%     FETCH by default is 100km.  U10 & FETCH must be equal sized vectors
%     or scalar, in which case S(:,i) gives the spectrum for U10(i) &
%     FETCH(i).
%
%    Notes:
%     - References:
%        +Pierson and Moskowitz, 1964, A Proposed Spectral Form for Fully
%         Developed Wind Seas Based on the Similarity Theory of S. A.
%         Kitaigorodskii, JGR, 69, 5181—5203
%        +Hasselmann, Barnett, Bouws, Carlson, Cartwright, Enke, Ewing,
%         Gienapp, Hasselmann, Kruseman, Meerburg, Müller, Olbers, Richter,
%         Sell, and Walden, 1973, Measurements of wind-wave growth and
%         swell decay during the Joint North Sea Project (JONSWAP),
%         Ergänzungsheft zur Deutschen Hydrographischen Zeitschrift,
%         No. 12, A8
%     - Wave amplitude:
%
%                     _f2
%        2           |
%       A  =  2 *    | S df
%                   _|
%                 f1
%
%     - Pierson-Moskowitz significant wave height:
%
%                       2
%       H    = 0.22 * U   / g
%        1/3           10
%
%     - JONSWAP significant wave height:
%
%                          2
%       H    = 1.67e-7 * U   * FETCH / g
%        1/3              10
%
%    Examples:
%     % The effect of windspeed on the Pierson-Moskowitz spectrum:
%     [s,f]=sea_spectrum(20:-2.5:10);
%     figure; plot(f,s);
%     xlabel('Freq (Hz)'); ylabel('PSD (m^2/Hz)');
%     xlim([0 .3]); ylim([0 110]);
%     legend(strcat(num2str((20:-2.5:10)'),'m/s'));
%     title('Pierson-Moskowitz spectrum vs. windspeed');
%
%     % The effect of fetch on the JONSWAP spectrum
%     [s,f]=sea_spectrum(7,[],[40 30 20 15 10],'j');
%     figure; plot(f,s);
%     xlabel('Freq (Hz)'); ylabel('PSD (m^2/Hz)');
%     xlim([0 .7]); ylim([0 .7]);
%     legend(strcat(num2str([40 30 20 15 10]'),'km'));
%     title('JONSWAP spectrum vs. fetch for 7m/s windspeed');
%
%     % Reproduce Figure 2 of Bromirski & Duennebier 2002:
%     [s,f]=sea_spectrum([5 7.5 10 15 20 30]);
%     figure; plot(f,20*log10(s),'k','linewidth',2);
%     xlabel('Frequency (Hz)'); ylabel('dB (m^2/Hz)');
%     text(0.26,-25,'5'); text(0.17,-7,'7.5');
%     text(0.125,5,'10'); text(0.08,22,'15');
%     text(0.06,34,'20'); text(0.06,56,'30');
%     setfonts(gca,'fontweight','bold');
%     xlim([0 .5]); ylim([-40 60]); grid on;
%
%    See also: SWELL_FORWARD, SWELL_VELOCITY, SWELL_BACKPROJ, LOCATE_STORMS

%     Version History:
%        Feb. 11, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2013 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% defaults
if(nargin<2 || isempty(freq)); freq=10.^(-2:.001:0)'; end
if(nargin<3 || isempty(fetch)); fetch=100; end
if(nargin<4 || isempty(type)); type='PM'; end

% check inputs
if(~isnumeric(u10) || ~isreal(u10) || any(u10<0) || ~isvector(u10))
    error('seizmo:sea_spectrum:badInput',...
        'U10 must be a positive real value in m/s!');
elseif(~isnumeric(freq) || ~isreal(freq) || any(freq<0) || ~isvector(freq))
    error('seizmo:sea_spectrum:badInput',...
        'FREQ must be positive real values in Hz!');
elseif(~isnumeric(fetch) || ~isreal(fetch) ...
        || any(fetch<0) || ~isvector(fetch))
    error('seizmo:sea_spectrum:badInput',...
        'FETCH must be positive real values in km!');
elseif(~isscalar(u10) && ~isscalar(fetch) && numel(u10)~=numel(fetch))
    error('seizmo:sea_spectrum:badInput',...
        'Non-scalar U10 & FETCH must match in size!');
elseif(~ischar(type) || ~isvector(type) || size(type,1)~=1)
    error('seizmo:sea_spectrum:badInput',...
        'TYPE must be a string!');
end

% gravity
g=9.81; % m/s/s

% fetch (km => m)
fetch=fetch*1e3;

% expand fetch/u10
if(isscalar(u10) && numel(fetch)>1); u10=u10(ones(size(fetch))); end
if(isscalar(fetch) && numel(u10)>1); fetch=fetch(ones(size(u10))); end

% preallocate
s=nan(numel(freq),numel(u10));

% loop over u10/fetch
for i=1:numel(u10)
    % constants
    switch lower(type)
        case {'pm' 'p' 'm' 'pierson' 'moskowitz' 'pierson-moskowitz'}
            fp=0.855*g/u10(i); % rad/s
            alpha=0.0081; % Phillip's constant
            gamma=1;
        case {'j' 'js' 'jonswap'}
            fp=22*(g^2/(u10(i)*fetch(i)))^(1/3); % rad/s
            alpha=0.076*(u10(i)^2/(fetch(i)*g))^.22;
            gamma=3.3;
        otherwise
            error('seizmo:sea_spectrum:badInput',...
                'Unknown sea spectrum: %s !',type);
    end
    fp=fp/(2*pi); % rad/s => Hz
    
    % spectrum
    sigma=0.07*(freq<=fp)+0.09*(freq>fp);
    s(:,i)=(alpha*g^2./(2*pi*freq).^5)...
        .*exp(-5/4*(fp./freq).^4)...
        .*gamma.^exp(-(freq-fp).^2./(2*sigma.^2*fp^2));
end

% needed for rad/s => Hz conversion
s=2*pi*s;

end
