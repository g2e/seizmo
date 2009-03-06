function [data,CONF]=filter_postprep(data,CONF)
% reads in data records and preps them for aligning

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
    CONF.SAMPLERATE=CONF.RATE1;
end

% ANALYSIS GROUND UNITS (DATA ASSUMED TO BE IN DISPLACEMENT BEFOREHAND!)
if(strcmpi(CONF.GUNITS,'velo'))
    data=dif(data);
elseif(strcmpi(CONF.GUNITS,'accel'))
    data=dif(dif(data));
end

end
