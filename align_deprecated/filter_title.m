function CONF=filter_title(CONF)
%FILTER_TITLE    Makes title for figures based on filter

if(strcmp(CONF.FILTTYPE,'low'))
    CONF.TITLE=['LOWPASS ' upper(CONF.FILTSTYLE) ' FILTERED AT ' num2str(CONF.PASSHIGH) ' Hz'];
elseif(strcmp(CONF.FILTTYPE,'high'))
    CONF.TITLE=['HIGHPASS ' upper(CONF.FILTSTYLE) ' FILTERED AT ' num2str(CONF.PASSLOW) ' Hz'];
elseif(strcmp(CONF.TYPE,'bandpass'))
    CONF.TITLE=['BANDPASS ' upper(CONF.FILTSTYLE) ' FILTERED BETWEEN ' num2str(CONF.PASSLOW) ' AND ' num2str(CONF.PASSHIGH) ' Hz'];
elseif(strcmp(CONF.TYPE,'notch'))
    CONF.TITLE=['NOTCH ' upper(CONF.FILTSTYLE) ' FILTERED BETWEEN ' num2str(CONF.PASSLOW) ' AND ' num2str(CONF.PASSHIGH) ' Hz'];
else
    error('unknown filter type')
end

end