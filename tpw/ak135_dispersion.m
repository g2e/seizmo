function [phvel]=ak135_dispersion(freq)
%AK135_DISPERSION    Returns fundamental mode AK135 dispersion
%
%    Usage:    phvel=ak135_dispersion(freq)
%
%    Description:
%     PHVEL=AK135_DISPERSION(FREQ) returns the Rayleigh wave phase
%     velocities PHVEL at the frequencies specified in FREQ.
%
%    Notes:
%     - At present this function is quite crappy.
%
%    Examples:
%     % Plot up a AK135 dispersion curve:
%     freq=0.005:0.001:0.1;
%     phvel=ak135_dispersion(freq);
%     plot(freq,phvel)
%
%    See also: PREM_DISPERSION, AK135

%     Version History:
%        Apr. 22, 2010 - initial version
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 09:15 GMT

% todo
% - get better dispersion values using Saito code
% - Love waves
% - group velocities
% - quality factor

% period vs ak135 phase velocity
ak135=[ % from Saito code
15.00000 3.397792
20.00000 3.578033
25.00000 3.731125
30.00000 3.833094
35.00000 3.897960
40.00000 3.941497
45.00000 3.973063
50.00000 3.997697
55.00000 4.018161
60.00000 4.036049
65.00000 4.052334
70.00000 4.067642
75.00000 4.082396
80.00000 4.096896
85.00000 4.111366
90.00000 4.125973
95.00000 4.140852
100.00000 4.156106
105.00000 4.171822
110.00000 4.188069
115.00000 4.204904
120.00000 4.222375
125.00000 4.240523
130.00000 4.259378
135.00000 4.278970
140.00000 4.299321
145.00000 4.320450
150.00000 4.342374
155.00000 4.365106
160.00000 4.388659
165.00000 4.413038
170.00000 4.438251
175.00000 4.464303
180.00000 4.491195
185.00000 4.518925
190.00000 4.547490
195.00000 4.576883
200.00000 4.607094
];

% linear interpolation
phvel=interp1(ak135(:,1),ak135(:,2),1./freq,[],'extrap');

end
