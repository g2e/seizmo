function [data,CONF]=filter_preprep(data,CONF)
% reads in data records and preps them for aligning

% READ IN DATA RECORDS
if(CONF.PREWINDOW); data=rpdw(data,CONF.PREWINREF1,CONF.PREWINOFF1,CONF.PREWINREF2,CONF.PREWINOFF2,CONF.PREWINFILL,CONF.PREWINFILLER); %%%%
else data=rdata(data); end

% DOUBLE PRECISION
if(CONF.DP)
    data=fixdelta(doubleit(data));
end

% PRE-FILTER DETREND
if(CONF.RDRIFT0==3); data=rdrift(data);
elseif(CONF.RDRIFT0==2); data=rslope(data);
elseif(CONF.RDRIFT0==1); data=rmean(data);
end

% PRE-FILTER DEAD TRACE REMOVAL
if(CONF.RDEAD0)
    data=rdead(data);
end

% PRE-FILTER TAPER
if(CONF.TAPER0)
    data=taper(data,CONF.TAPERHW0,CONF.TAPERTYPE0,CONF.TAPEROPT0);
end

% PRE-FILTER RESAMPLE
if(CONF.RESAMPLE0)
    if(CONF.FIRINTERP0); data=syncsr(data,CONF.RATE0);
    else data=intrpol8(data,CONF.RATE0,CONF.INTERPMETH0); end
    CONF.SAMPLERATE=CONF.RATE0;
end

end
