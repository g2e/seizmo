% How to create the IRIS SAC PoleZero Database:
% 1. Create directory of all sac polezero files
%   a. retreive all dataless (/sks/ggeuler/CATALOGS/scripts/cron/dataless)
%   b. extract all sac polezeros (/sks/ggeuler/CATALOGS/make_sacpz)
% 2. Build in matlab:
%        net=xdir('/sks/ggeuler/CATALOGS/sacpz');
%        net(strcmp({net.name},'.') | strcmp({net.name},'..'))=[];
%        matlabpool(8);
%        sacpzdb_tmp=cell(numel(net),1);
%        name=cell(numel(net),1);
%        parfor i=1:numel(net)
%            disp(net(i).name);
%            tmp=getwords(net(i).name,'.'); % remove years off of temp array directories
%            name{i}=tmp{1};                % to get the network name (required)
%            sacpzdb_tmp{i}=makesacpzdb(fullfile(net(i).path,net(i).name));
%        end
%        sacpzdb_tmp=cat(1,sacpzdb_tmp{:});
%        matlabpool close;
% 3. Clear the path: [sacpzdb_tmp.path]=deal('./');
% 4. In matlab (breaking up db into sub-dbs by network):
%    knetwk={sacpzdb_tmp.knetwk};
%    nets=unique(knetwk);
%    for i=1:numel(nets)
%        % handle networks that begin with a digit
%        if(isstrprop(nets{i}(1),'digit'))
%            netname=['A_' nets{i}];
%        else
%            netname=nets{i};
%        end
%        sacpzdb.(netname)=sacpzdb_tmp(strcmpi(knetwk,nets(i)));
%    end
% 5. Run cplxpair fixes in matlab: sacpzdb=fix_bad_sacpz(sacpzdb);
% 6. Save database in matlab: save sacpzdb -struct sacpzdb
% 7. Move sacpzdb.mat to seizmo/response directory

% Last Updated: March 5, 2011 by Garrett Euler
