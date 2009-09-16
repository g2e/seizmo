function [head]=vf_ch_nzmonth(def,head,value)
%VF_CH_NZMONTH    Sets virtual field NZMONTH

% if month is undef => set jday undef
% if year is undef => force year to be undef on return
% no matter month => set cday to 1

value(:,2)=1; % set calendar day to the first of the month
head(def.reftime(1:2),:)=cal2doy([head(def.reftime(1),:).' value]).';

end
