function []=iris_sacpzdb_build(path)
%IRIS_SACPZDB_BUILD    Creates the IRIS SAC PoleZero Database
%
%    Usage:    iris_sacpzdb_build
%              iris_sacpzdb_build(path)
%
%    Description:
%     IRIS_SACPZDB_BUILD creates the polezero database that contains the
%     polezero responses of most seismometers with data available from
%     IRIS (Incorporated Research Institutions for Seismology).  This db is
%     used by SEIZMO to automatically determine the instrument response of
%     records.  The db is more complicated than what MAKESACPZDB creates
%     (which is probably what you should be looking at).  Anyway, this
%     assumes a specific directory structure that is made by some shell
%     scripts I wrote that are not included in SEIZMO so you need those for
%     this to do anything useful.
%
%     IRIS_SACPZDB_BUILD(PATH) allows changing the path to the specific
%     directory structure.  The default assumes a location on a WashU
%     seismology group server.
%
%    Notes:
%     - Talk to me before you use this.
%     - Uses parallel computing toolbox.
%
%    Examples:
%     % This is what I do from a shell:
%     /sks/ggeuler/CATALOGS/scripts/cron/dataless
%     /sks/ggeuler/CATALOGS/make_sacpz
%     matlabcli -r 'iris_sacpzdb_build;exit'
%     cp sacpzdb_fixed.mat /opt/seizmo/seizmo/response/sacpzdb.mat
%
%    See also: MAKESACPZDB, FIX_BAD_SACPZ

%     Version History:
%        Dec.  2, 2009 - initial version
%        May  28, 2010 - handle networks starting with a digit
%        Mar.  5, 2011 - added parallelization
%        Dec. 22, 2011 - convert iris_sacpzdb_howto notes to this script
%        Feb.  3, 2012 - made documentation, removed inert code, use new
%                        form of bad_sacpz_cmplx (avoids crashes)
%        Mar. 15, 2012 - fix parallel verbose
%        Mar. 28, 2013 - turn verbosity back, slight change in output names
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 28, 2013 at 11:15 GMT

% todo:

% default path
if(nargin<1 || isempty(path)); path='/sks/ggeuler/CATALOGS/sacpz'; end

% build db in matlab
net=xdir(path);
net(strcmp({net.name},'.') | strcmp({net.name},'..'))=[];
verbose=seizmoverbose;
matlabpool(8);
sacpzdb_tmp=cell(numel(net),1);
parfor i=1:numel(net)
    disp(net(i).name);
    seizmoverbose(false);
    sacpzdb_tmp{i}=makesacpzdb(fullfile(net(i).path,net(i).name));
end
matlabpool close;
seizmoverbose(verbose);

% flatten db & clear paths
disp('FLATTENING DB');
sacpzdb_tmp=cat(1,sacpzdb_tmp{:});
[sacpzdb_tmp.path]=deal('./');

% break into subdbs (for speed)
disp('SPLITTING INTO SUBDBs');
knetwk={sacpzdb_tmp.knetwk};
nets=unique(knetwk);
for i=1:numel(nets)
    % handle networks that begin with a digit
    if(isstrprop(nets{i}(1),'digit'))
        netname=['A_' nets{i}];
    else
        netname=nets{i};
    end
    sacpzdb.(netname)=sacpzdb_tmp(strcmpi(knetwk,nets(i)));
end

% save unfixed db (just in case)
disp('FIXING KNOWN BAD RESPONSES & SAVING');
save sacpzdb_unfixed -struct sacpzdb

% apply fixes and save
sacpzdb=fix_bad_sacpz(sacpzdb);
[badpz,sacpzdb]=bad_sacpz_cplxpair(sacpzdb);
save sacpzdb_fixed -struct sacpzdb

end
