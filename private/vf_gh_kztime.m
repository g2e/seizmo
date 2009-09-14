function [value]=vf_gh_kztime(def,head)
%VF_GH_KZTIME    Returns value for virtual field KZTIME

value=head(def.reftime(3:6),:);
value=sprintf('%02d:%02d:%02d.%03d',value);
value=cellstr(reshape(value,12,[]).');

end
