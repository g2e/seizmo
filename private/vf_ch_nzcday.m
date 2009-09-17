function [head]=vf_ch_nzcday(def,head,value)
%VF_CH_NZCDAY    Sets virtual field NZCDAY

% pull in current year & jday
yj=head(def.reftime(1:2),:).';

% who's (un)defined
bady=yj(:,1)==def.undef.ntype;
badj=yj(:,2)==def.undef.ntype;
badc=value==def.undef.ntype;
goodc=~badc;

% get current year
dt=datevec(now);

% undef jday if cday undef
if(any(~badj & badc))
    yj(~badj & badc,2)=def.undef.ntype;
end

% set undef year to current if cday ok
if(any(bady & goodc))
    yj(bady & goodc,1)=dt(1);
end

% set undef jday to 1 if cday ok
if(any(badj & goodc))
    yj(badj & goodc,2)=1;
end

% get new jday based on cday
if(any(goodc))
    cal=doy2cal(yj(goodc,:));
    cal(:,3)=value(goodc);
    yj(goodc,:)=cal2doy(cal);
end

% set header
head(def.reftime(1:2),:)=yj.';

end
