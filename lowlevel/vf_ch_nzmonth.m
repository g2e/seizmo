function [head]=vf_ch_nzmonth(def,head,value)
%VF_CH_NZMONTH    Sets virtual field NZMONTH

% pull in current year & jday
yj=head(def.reftime(1:2),:).';

% who's (un)defined
bady=yj(:,1)==def.undef.ntype;
badj=yj(:,2)==def.undef.ntype;
badm=value==def.undef.ntype;
goodm=~badm;

% get current year
dt=datevec(now);

% undef jday if month undef
if(any(~badj & badm))
    yj(~badj & badm,2)=def.undef.ntype;
end

% set undef year to current if month ok
if(any(bady & goodm))
    yj(bady & goodm,1)=dt(1);
end

% get new jday based on month
if(any(goodm))
    cal=[yj(goodm,1) value(goodm)];
    cal(:,3)=1; % cday always set to 1 if month set
    yj(goodm,:)=cal2doy(cal);
end

% set header
head(def.reftime(1:2),:)=yj.';

end
