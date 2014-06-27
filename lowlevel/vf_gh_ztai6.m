function [value,good]=vf_gh_ztai6(def,head)
%VF_GH_ZTAI6    Returns value for virtual field ZTAI6

% get reference time
tmp=head(def.reftime,:);

% who's  (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [false(1,nv); (tmp(2,:)<1 | tmp(2,:)>366); ...
    (tmp(3,:)<0 | tmp(3,:)>23); (tmp(4,:)<0 | tmp(4,:)>59); ...
    (tmp(5,:)<0 | tmp(5,:)>60); (tmp(6,:)<0 | tmp(6,:)>1000)])==0;

% default [yr jday hr mn secs] all undef
value(nv,6)=def.undef.ntype;

if(any(good))
    % hr, min already known
    value(good,4:5)=tmp(3:4,good).';
    
    % get month/cday
    value(good,1:3)=doy2cal(tmp(1:2,good).');
    
    % get secs
    value(good,6)=(tmp(5,good)+tmp(6,good)/1000).';
    
    % convert to tai
    value(good,:)=utc2tai(value(good,:));
end

% normally wrap in cell
if(nargout==1); value=mat2cell(value,ones(nv,1)); end

end
