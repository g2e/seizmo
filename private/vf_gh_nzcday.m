function [value]=vf_gh_nzcday(def,head)
%VF_GH_NZCDAY    Returns value for virtual field NZCDAY

value=doy2cal(head(def.reftime(1:2),:).');
value=value(:,3);

end
