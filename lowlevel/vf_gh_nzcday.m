function [value]=vf_gh_nzcday(def,head)
%VF_GH_NZCDAY    Returns value for virtual field NZCDAY

% get current reftime
tmp=head(def.reftime(1:2),:);

% who's (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [false(1,nv); (tmp(2,:)<1 | tmp(2,:)>366)])==0;

% default to all undef
value(nv,1)=def.undef.ntype;

if(any(good))
    % get month/cday
    tmp=doy2cal(tmp(:,good).');
    
    % separate and add into output
    value(good,1)=tmp(:,3);
end

end
