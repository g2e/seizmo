function [value]=vf_gh_kzdate(def,head)
%VF_GH_KZDATE    Returns value for virtual field KZDATE

value=head(def.reftime(1:2),:);
cal=doy2cal(value(1:2,:).');
value=sprintf('%04d-%02d-%02d (%03d)',[cal.'; value(2,:)]);
value=cellstr(reshape(value,16,[]).');

end
