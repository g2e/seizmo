function [head]=vf_ch_nzcday(def,head,value)
%VF_CH_NZCDAY    Sets virtual field NZCDAY

yj=head(def.reftime(1:2),:).';
bady=yj(:,1)==def.undef.ntype;
badj=yj(:,2)==def.undef.ntype;
bad2=bady & badj;
dt=datevec(now); dt2=cal2doy(dt(1:3));
if(any(bad2))
    % both year and jday undef
    % - use current year/doy
    yj(bad2,1)=dt(1);
    yj(bad2,2)=dt2(2);
end
if(any(bady & ~badj))
    % just year undef
    % - use current year
    yj(bady & ~badj,1)=dt(1);
end
if(any(badj & ~bady))
    % just jday undef
    % - set jday as 1
    yj(badj & ~bady,2)=1;
end

cal=doy2cal(yj);
cal(:,3)=value;
head(def.reftime(1:2),:)=cal2doy(cal).';

end
