function [value]=vf_lh_ztai(def,head)
%VF_LH_ZTAI    Returns value for virtual field ZTAI as a string

% get current reftime
tmp=head(def.reftime,:);

% who's (un)defined
nv=size(head,2);
good=sum(isnan(tmp) | isinf(tmp) | tmp==def.undef.ntype ...
    | tmp~=round(tmp) | [false(1,nv); (tmp(2,:)<1 | tmp(2,:)>366); ...
    (tmp(3,:)<0 | tmp(3,:)>23); (tmp(4,:)<0 | tmp(4,:)>59); ...
    (tmp(5,:)<0 | tmp(5,:)>60); (tmp(6,:)<0 | tmp(6,:)>1000)])==0;

% default to undef
value(nv,1)={def.undef.stype};

if(any(good))
    % get secs
    tmp(7,good)=(tmp(5,good)+tmp(6,good)/1000).';
    
    % convert to tai
    tmp([1:4 7],good)=utc2tai(tmp([1:4 7],good).').';
    
    % get sec/msec
    tmp(7,good)=round(1000*tmp(7,good));
    tmp(5,good)=fix(tmp(7,good)/1000);
    tmp(6,good)=mod(tmp(7,good),1000);
    
    % get month/cday
    cal=doy2cal(tmp(1:2,good).');
    
    % make string
    tmp=sprintf('%04d-%02d-%02d (%03d) %02d:%02d:%02d.%03d',...
        [cal.'; tmp(2:6,good)]);
    
    % separate and add into output
    value(good,1)=cellstr(reshape(tmp,29,[]).');
end

end
