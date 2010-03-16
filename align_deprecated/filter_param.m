function CONF=filter_param(CONF)
%FILTER_PARAM    Configures some parameters automatically based on filter

% input for iirfilter
CONF.FILTATTEN=[CONF.FILTSTOPATTEN CONF.FILTPASSATTEN];
if(strcmp(CONF.FILTTYPE,'low'))
    CONF.FILTLIMITS=[CONF.PASSHIGH CONF.STOPHIGH];
elseif(strcmp(CONF.FILTTYPE,'high'))
    CONF.FILTLIMITS=[CONF.STOPLOW CONF.PASSLOW];
elseif(strcmp(CONF.FILTTYPE,'bandpass'))
    CONF.FILTLIMITS=[CONF.STOPLOW CONF.PASSLOW CONF.PASSHIGH CONF.STOPHIGH];
elseif(strcmp(CONF.FILTTYPE,'notch'))
    CONF.FILTLIMITS=[CONF.STOPLOW CONF.PASSLOW CONF.PASSHIGH CONF.STOPHIGH];
else
    error('unknown filter type')
end

% input for xc multi-picker
if(isempty(CONF.SPACING)); CONF.SPACING=ceil(CONF.SAMPLERATE/(4*CONF.PASSHIGH)); end
if(isempty(CONF.ADJACENT)); CONF.ADJACENT=floor(CONF.SAMPLERATE/(16*CONF.PASSHIGH)); end

end