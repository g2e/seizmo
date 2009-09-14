function [head]=vf_ch_nzcday(def,head,value)
%VF_CH_NZCDAY    Sets virtual field NZCDAY

cal=doy2cal(head(def.reftime(1:2),:).');
cal(:,3)=value;
head(def.reftime(1:2),:)=cal2doy(cal).';

end
