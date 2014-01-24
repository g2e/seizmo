function [mod]=readcrust2(keyfile,typefile,elevfile)
%READCRUST2    Reads Crust2.0 files, putting them into a struct
%
%    Usage:    mod=readcrust2(key,type,elev)
%
%    Description:
%     MOD=READCRUST2(KEY,TYPE,ELEV) reads in the Crust2.0 model from the
%     files given by KEY, TYPE & ELEV.  Note that the output is a struct
%     with the model values stored in integer format for memory savings.
%     This format needs to be converted to double precision and divided by
%     100 to get the actual values.  Please use GETCRUST to access the
%     values as this is not really meant to be used for routine access.
%
%    Notes:
%     - Also works for Crust5.1
%
%    Examples:
%     % To read in Crust2.0 from an expanded archive:
%     mod=readcrust2('CNtype2_key.txt','CNtype2.txt','CNelevatio2.txt');
%
%     % To read in Crust5.1:
%     mod=readcrust2('CNtype_key.txt','CNtype.txt','CNelevatio.txt');
%
%    See also: READCRUST10, GETCRUST

%     Version History:
%        May  18, 2010 - initial version
%        Aug.  5, 2010 - minor doc update
%        Jan. 23, 2014 - added full docs, changed to expanded output rather
%                        than type/key for efficiency with GETCRUST update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 02:45 GMT

% todo:

% read in key file
lines=getwords(readtxt(keyfile),sprintf('\n'));
lines=lines(6:end); % skip header
ntypes=numel(lines)/5;
print_time_left(0,ntypes+7);
for i=1:ntypes
    words=getwords(lines{(i-1)*5+1});
    key(i).type=words{1};
    key(i).desc=joinwords(words(2:end));
    key(i).vp=str2double(getwords(lines{(i-1)*5+2}));
    key(i).vs=str2double(getwords(lines{(i-1)*5+3}));
    key(i).rho=str2double(getwords(lines{(i-1)*5+4}));
    key(i).thick=str2double(getwords(lines{(i-1)*5+5}));
    key(i).thick=key(i).thick(1:end-2);
    print_time_left(i,ntypes+7);
end
tk.key=key;

% read in topo grid
lines=getwords(readtxt(elevfile),sprintf('\n'));
tk.elev=str2num(char(lines(2:end))); % drops lons
tk.elev=tk.elev(:,2:end); % drops lats
tk.elev=round(tk.elev/10)*10; % remove 1 sig. digit (elev in 10m steps)
print_time_left(ntypes+1,ntypes+7);

% read in type grid
sz=size(tk.elev);
lines=getwords(readtxt(typefile),sprintf('\n'));
type=reshape(getwords(joinwords(lines(2:end),' ')),[sz(2)+1 sz(1)])';
type=type(:,2:end); % drop lats
[idx,idx]=ismember(type,{tk.key.type}');
tk.type=idx;
print_time_left(ntypes+2,ntypes+7);

% add in moho info
%tk.moho=nan(sz);
%for i=1:prod(sz)
%    tk.moho(i)=tk.elev(i)/1000 ...
%        -sum(tk.key(tk.type(i)).thick([1 3:end]));
%end
%tk.moho=-1*tk.moho;
%print_time_left(ntypes+3,ntypes+3);

% convert to expanded form
tmp=tk.key(tk.type(:));
mod.top=tk.elev(:)/1000;
mod.vp=cat(1,tmp.vp);
mod.vs=cat(1,tmp.vs);
mod.rho=cat(1,tmp.rho);
mod.thk=cat(1,tmp.thick);
print_time_left(ntypes+3,ntypes+7);

% switch water/ice & add middle seds layer
mod.vp=mod.vp(:,[2 1 3 3 4:8]);
mod.vp(:,4)=0;
mod.vs=mod.vs(:,[2 1 3 3 4:8]);
mod.vs(:,4)=0;
mod.rho=mod.rho(:,[2 1 3 3 4:8]);
mod.rho(:,4)=0;
mod.thk=mod.thk(:,[2 1 3 3 4:7]);
mod.thk(:,4)=0;
print_time_left(ntypes+4,ntypes+7);

% use elevation and thickness to get top
mod.top=mod.top(:,ones(1,9));
mod.top(:,3:9)=mod.top(:,3:9)-cumsum(mod.thk(:,2:8),2);
mod.top(mod.top(:,1)<0,1)=0; % death valley problem...
print_time_left(ntypes+5,ntypes+7);

% remove thickness
mod=rmfield(mod,'thk');
print_time_left(ntypes+6,ntypes+7);

% convert to integers for efficient storage
mod.top=int16(mod.top*100);
mod.rho=int16(mod.rho*100);
mod.vp=int16(mod.vp*100);
mod.vs=int16(mod.vs*100);
print_time_left(ntypes+7,ntypes+7);

end
