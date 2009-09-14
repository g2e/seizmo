function [head]=vf_ch_nzmonth(def,head,value)
%VF_CH_NZMONTH    Sets virtual field NZMONTH

value(:,2)=1; % set calendar day to the first of the month
head(def.reftime(1:2),:)=cal2doy([head(def.reftime(1),:).' value]).';

end
