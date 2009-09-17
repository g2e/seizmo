function [head]=vf_ch_kzdttm(def,head,value)
%VF_CH_KZDTTM    Sets virtual field KZDTTM

% number of records
nv=numel(value);

% default output
dttm=def.undef.ntype*ones(nv,6);
good=false(nv,1);

% who's (un)defined
num=true(1,29);
num(1,[5 8 11 12 16 17 20 23 26])=false;
for i=1:nv
    % must be a 29 character string with proper separators
    % and numeric chars elsewhere
    if(ischar(value{i}) && isequal(size(value{i}),[1 29]) ...
            && strcmp(value{i}(~num),'-- () ::.') ...
            && all(isstrprop(value{i}(num),'digit')))
        good(i)=true;
    end
end

% parse good
if(any(good))
    value=char(value(good,1));
    value=[cellstr(value(:,1:4)) cellstr(value(:,13:15)) ...
        cellstr(value(:,18:19)) cellstr(value(:,21:22)) ...
        cellstr(value(:,24:25)) cellstr(value(:,27:29))];
    dttm(good,:)=str2double(value);
end

% set header
head(def.reftime,:)=dttm.';

end
