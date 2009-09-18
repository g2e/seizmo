function [head]=vf_ch_z(def,head,value)
%VF_CH_Z    Sets virtual field Z

% number of records
nv=numel(value);

% default to undef
head(def.reftime,:)=def.undef.ntype;
good5=false(nv,1);
good6=good5;

% who's (un)defined
for i=1:nv
    % must be 1x5 or 1x6 numeric array
    if(isnumeric(value{i}))
        if(isequal(size(value{i}),[1 5]) ...
                && ~any(value{i}==def.undef.ntype) ...
                && ~any(isnan(value{i}) | isinf(value{i})) ...
                && isequal(value{i}(1:4),round(value{i}(1:4))))
            good5(i)=true;
        elseif(isequal(size(value{i}),[1 6]) ...
                && ~any(value{i}==def.undef.ntype) ...
                && ~any(isnan(value{i}) | isinf(value{i})) ...
                && isequal(value{i}(1:5),round(value{i}(1:5))))
            good6(i)=true;
        end
    end
end

% skip empty
if(any(good5))
    % fix times
    value5=fixtimes(cell2mat(value(good5,1)));
    
    % get sec and msec
    value5(:,5)=round(1000*value5(:,5));
    value5(:,6)=mod(value5(:,5),1000);
    value5(:,5)=fix(value5(:,5)/1000);
    
    % set header
    head(def.reftime,good5)=value5.';
end
if(any(good6))
    % fix times
    value6=fixtimes(cell2mat(value(good6,1)));
    
    % convert month/cday to jday
    value=cal2doy(value6(:,1:3));
    value6=[value value6(:,4:6)];
    
    % get sec and msec
    value6(:,5)=round(1000*value6(:,5));
    value6(:,6)=mod(value6(:,5),1000);
    value6(:,5)=fix(value6(:,5)/1000);
    
    % set header
    head(def.reftime,good6)=value6.';
end

end
