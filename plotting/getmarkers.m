function [names,times]=getmarkers(data)
%GETMARKERS    Gets marker info for SEIZMO plotting

% get marker name/time
[ka,a,kf,f,ko,o,kt,t]=getheader(data,'ka','a','kf','f','ko','o','kt','t');

% fix names
badka=strcmpi(ka,'nan');
if(any(badka)); ka(badka)={'a'}; end
badkf=strcmpi(kf,'nan');
if(any(badkf)); kf(badkf)={'f'}; end
badko=strcmpi(ko,'nan');
if(any(badko)); ko(badko)={'o'}; end
badkt=strcmpi(kt,'nan');
if(any(badkt))
    [row,col]=find(badkt);
    kt(badkt)=strcat('t',cellstr(num2str(row)));
end

% combine under single arrays
names=[ka kf ko kt];
times=[a f o t];

end
