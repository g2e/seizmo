function [head]=vf_ch_kztime(def,head,value)
%VF_CH_KZTIME    Sets virtual field KZTIME

% number of records
nv=numel(value);

% default output
tm=def.undef.ntype*ones(nv,4);
good=false(nv,1);

% who's (un)defined
num=true(1,12);
num(1,[3 6 9])=false;
for i=1:nv
    % must be a 12 character string with proper separators
    % and numeric chars elsewhere
    if(ischar(value{i}) && isequal(size(value{i}),[1 12]) ...
            && strcmp(value{i}(~num),'::.') ...
            && all(isstrprop(value{i}(num),'digit')))
        good(i)=true;
    end
end

% parse good
if(any(good))
    value=char(value(good,1));
    value=[cellstr(value(:,1:2)) cellstr(value(:,4:5)) ...
        cellstr(value(:,7:8)) cellstr(value(:,10:12))];
    tm(good,:)=str2double(value);
end

% set header
head(def.reftime(3:6),:)=tm.';

end
