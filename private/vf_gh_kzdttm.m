function [value]=vf_gh_kzdttm(def,head)
%VF_GH_KZDTTM    Returns value for virtual field KZDTTM

value=head(def.reftime(1:6),:);
cal=doy2cal(value(1:2,:).');
value=sprintf('%04d-%02d-%02d (%03d) %02d:%02d:%02d.%03d',...
    [cal.'; value(2:6,:)]);
value=cellstr(reshape(value,29,[]).');

end
