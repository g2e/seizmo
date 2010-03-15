function [value]=vf_gh_kzdttm(def,head)
%VF_GH_KZDTTM    Returns value for virtual field KZDTTM

% get current reftime
tmp=head(def.reftime(1:6),:);

% who's (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [false(1,nv); (tmp(2,:)<1 | tmp(2,:)>366); ...
    (tmp(3,:)<0 | tmp(3,:)>23); (tmp(4,:)<0 | tmp(4,:)>59); ...
    (tmp(5,:)<0 | tmp(5,:)>60); (tmp(6,:)<0 | tmp(6,:)>1000)])==0;

% default to all undef
value(nv,1)={def.undef.stype};

if(any(good))
    % get month/cday
    cal=doy2cal(tmp(1:2,good).');
    
    % make string
    tmp=sprintf('%04d-%02d-%02d (%03d) %02d:%02d:%02d.%03d',...
        [cal.'; tmp(2:6,good)]);
    
    % separate and add into output
    value(good,1)=cellstr(reshape(tmp,29,[]).');
end

end
