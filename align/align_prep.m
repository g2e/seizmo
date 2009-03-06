function [data]=align_prep(data,CONF)
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
    if(CONF.INTERP0); data=intrpol8(data,CONF.RATE0,CONF.INTERPMETH0);
    else data=syncsr(data,CONF.RATE0); end
end

% FILTER
if(CONF.FILTER)
    data=iirfilter(data,CONF.TYPE,CONF.STYLE,CONF.LIMITS,CONF.ORDER,CONF.NPASS,CONF.RIPPLE);
end


% POST-FILTER TAPER1
if(CONF.TAPER1)
    data=taper(data,CONF.TAPERHW1,CONF.TAPERTYPE1,CONF.TAPEROPT1);
end

% POST-FILTER DETREND
if(CONF.RDRIFT1==3); data=rdrift(data);
elseif(CONF.RDRIFT1==2); data=rslope(data);
elseif(CONF.RDRIFT1==1); data=rmean(data);
end

% POST-FILTER TAPER2
if(CONF.TAPER2)
    data=taper(data,CONF.TAPERHW2,CONF.TAPERTYPE2,CONF.TAPEROPT2);
end

% PRE-FILTER DEAD TRACE REMOVAL
if(CONF.RDEAD1)
    data=rdead(data);
end

% POST-FILTER RESAMPLE
if(CONF.RESAMPLE1)
    if(CONF.INTERP1); data=intrpol8(data,CONF.RATE1,CONF.INTERPMETH1);
    else data=syncsr(data,CONF.RATE1); end
end

% ANALYSIS GROUND UNITS (DATA EXPECTED TO BE IN DISPLACEMENT BEFOREHAND!)
if(strcmpi(CONF.GUNITS,'velo'))
    data=dif(data);
elseif(strcmpi(CONF.GUNITS,'accel'))
    data=dif(dif(data));
end

% ANALYZE ENVELOPES?
if(CONF.ENV)
    data=envelope(data);
end

end
