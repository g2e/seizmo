function [mod]=readcrust2(keyfile,typefile,elevfile)
%READCRUST2    Reads Crust2.0 files, putting them into a struct
%
% Usage: mod=readcrust2('CNtype2_key.txt','CNtype2.txt','CNelevatio2.txt')

% read in key file
lines=getwords(readtxt(keyfile),sprintf('\n'));
lines=lines(6:end); % skip header
ntypes=numel(lines)/5;
print_time_left(0,ntypes+2);
for i=1:ntypes
    words=getwords(lines{(i-1)*5+1});
    key(i).type=words{1};
    key(i).desc=joinwords(words(2:end));
    key(i).vp=str2double(getwords(lines{(i-1)*5+2}));
    key(i).vs=str2double(getwords(lines{(i-1)*5+3}));
    key(i).rho=str2double(getwords(lines{(i-1)*5+4}));
    key(i).thick=str2double(getwords(lines{(i-1)*5+5}));
    key(i).thick=key(i).thick(1:end-2);
    print_time_left(i,ntypes+3);
end
mod.key=key;

% read in topo grid
lines=getwords(readtxt(elevfile),sprintf('\n'));
mod.elev=str2num(char(lines(2:end)));
mod.elev=mod.elev(:,2:end); % drop lats
print_time_left(ntypes+1,ntypes+3);

% read in type grid
sz=size(mod.elev);
lines=getwords(readtxt(typefile),sprintf('\n'));
type=reshape(getwords(joinwords(lines(2:end),' ')),[sz(2)+1 sz(1)])';
type=type(:,2:end); % drop lats
[idx,idx]=ismember(type,{mod.key.type}');
mod.type=idx;
print_time_left(ntypes+2,ntypes+3);

% add in moho info
mod.moho=nan(sz);
for i=1:prod(sz)
    mod.moho(i)=mod.elev(i)/1000 ...
        -sum(mod.key(mod.type(i)).thick([1 3:end]));
end
mod.moho=-1*mod.moho;
print_time_left(ntypes+3,ntypes+3);

end
