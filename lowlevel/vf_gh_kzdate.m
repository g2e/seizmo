function [value]=vf_gh_kzdate(def,head)
%VF_GH_KZDATE    Returns value for virtual field KZDATE

% get current reftime
tmp=head(def.reftime(1:2),:);

% who's (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [false(1,nv); (tmp(2,:)<1 | tmp(2,:)>366)])==0;

% default to all undef
value(nv,1)={def.undef.stype};

if(any(good))
    % get month/cday
    cal=doy2cal(tmp(:,good).');
    
    % make string
    tmp=sprintf('%04d-%02d-%02d (%03d)',[cal.'; tmp(2,good)]);
    
    % separate and add into output
    value(good,1)=cellstr(reshape(tmp,16,[]).');
end

end
