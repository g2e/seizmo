% How to create the IRIS SAC PoleZero Database:
% 1. Create directory of all sac polezero files
%   a. retreive all dataless (/research/scripts/cron/dataless)
%   b. extract all sac polezeros (/research/catalogs/make_sacpz.bash)
% 2. Build in matlab: sacpzdb_tmp=makesacpzdb('/research/catalogs/sacpz');
% 3. In matlab (breaking up db into sub-dbs by network):
%    knetwk={sacpzdb_tmp.knetwk};
%    nets=unique(knetwk);
%    for i=1:numel(nets)
%        sacpzdb.(nets{i})=sacpzdb_tmp(strcmpi(knetwk,nets(i)));
%    end
% 4. Save database in Matlab: save sacpzdb -struct sacpzdb
% 5. Move sacpzdb.mat to seizmo/response directory