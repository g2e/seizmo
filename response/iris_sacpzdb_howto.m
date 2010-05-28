% How to create the IRIS SAC PoleZero Database:
% 1. Create directory of all sac polezero files
%   a. retreive all dataless (/p/ggeuler/scripts/cron/dataless)
%   b. extract all sac polezeros (/p/ggeuler/catalogs/make_sacpz.bash)
% 2. Build in matlab: sacpzdb_tmp=makesacpzdb('/p/ggeuler/catalogs/sacpz');
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

% Last Updated: May 28, 2010 by Garrett Euler