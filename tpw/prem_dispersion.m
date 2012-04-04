function [phvel]=prem_dispersion(freq)
%PREM_DISPERSION    Returns fundamental mode PREM dispersion
%
%    Usage:    phvel=prem_dispersion(freq)
%
%    Description:
%     PHVEL=PREM_DISPERSION(FREQ) returns the Rayleigh wave phase
%     velocities PHVEL at the frequencies specified in FREQ.
%
%    Notes:
%     - At present this function is crap.
%
%    Examples:
%     % Plot up a PREM dispersion curve:
%     freq=0.005:0.001:0.1;
%     phvel=prem_dispersion(freq);
%     plot(freq,phvel)
%
%    See also: AK135_DISPERSION, PREM

%     Version History:
%        Feb.  5, 2010 - initial version
%        Apr. 22, 2010 - slightly less crappy (use output from Saito Code)
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 09:15 GMT

% todo
% - get better dispersion values using Saito code
% - Love waves
% - group velocities
% - quality factor

% period vs prem phase velocity
prem=[ % from Saito code
  15.  3.48093891
  20.  3.79185605
  25.  3.89421678
  30.  3.94230103
  35.  3.97109962
  40.  3.9914825
  45.  4.00779486
  50.  4.02212572
  55.  4.03561306
  60.  4.04892397
  65.  4.06245136
  70.  4.07642984
  75.  4.09099054
  80.  4.1061964
  85.  4.12207556
  90.  4.13863087
  95.  4.15585232
  100.  4.17372847
  105.  4.19224453
  110.  4.21138763
  115.  4.23114777
  120.  4.25151777
  125.  4.27249527
  130.  4.29407978
  135.  4.31627369
  140.  4.33908224
  145.  4.36250973
  150.  4.38656569
  155.  4.41125822
  160.  4.43659496
  165.  4.46258497
  170.  4.48923683
  175.  4.51655626
  180.  4.54454947
  185.  4.57322025
  190.  4.60256863
  195.  4.63259411
  200.  4.66329098
];

% linear interpolation
phvel=interp1(prem(:,1),prem(:,2),1./freq,[],'extrap');

end
