function [value]=vf_gh_nzmonth(def,head)
%VF_GH_NZMONTH    Returns value for virtual field NZMONTH

value=doy2cal(head(def.reftime(1:2),:).');
value=value(:,2);

end
