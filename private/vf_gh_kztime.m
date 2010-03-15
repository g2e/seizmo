function [value]=vf_gh_kztime(def,head)
%VF_GH_KZTIME    Returns value for virtual field KZTIME

% get current reftime
tmp=head(def.reftime(3:6),:);

% who's (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [(tmp(1,:)<0 | tmp(1,:)>23); ...
    (tmp(2,:)<0 | tmp(2,:)>59); (tmp(3,:)<0 | tmp(3,:)>60); ...
    (tmp(4,:)<0 | tmp(4,:)>1000)])==0;

% default to all undef
value(nv,1)={def.undef.stype};

if(any(good))
    % make string
    tmp=sprintf('%02d:%02d:%02d.%03d',tmp(:,good));
    
    % separate and add into output
    value(good,1)=cellstr(reshape(tmp,12,[]).');
end

end
