function [head]=vf_ch_kzdate(def,head,value)
%VF_CH_KZDATE    Sets virtual field KZDATE

% number of records
nv=numel(value);

% default output
dt=def.undef.ntype*ones(nv,2);
good=false(nv,1);

% who's (un)defined
num=true(1,16);
num(1,[5 8 11 12 16])=false;
for i=1:nv
    % must be a 16 character string with proper separators
    % and numeric chars elsewhere
    if(ischar(value{i}) && isequal(size(value{i}),[1 16]) ...
            && strcmp(value{i}(~num),'-- ()') ...
            && all(isstrprop(value{i}(num),'digit')))
        good(i)=true;
    end
end

% parse good
if(any(good))
    value=char(value(good,1));
    value=[cellstr(value(:,1:4)) cellstr(value(:,13:15))];
    dt(good,:)=str2double(value);
end

% set header
head(def.reftime(1:2),:)=dt.';

end
