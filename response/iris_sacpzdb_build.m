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
%     /opt/seizmo/sacpzdb/cron/dataless
%     /opt/seizmo/sacpzdb/make_sacpz
%     matlabcli -r 'iris_sacpzdb_build;exit'
%     cp sacpzdb.mat /opt/seizmo/seizmo/response/
%     cd /opt/seizmo/seizmo
%     zip seizmo_iris_sacpzdb.zip response/sacpzdb.mat
%
%    See also: MAKESACPZDB, FIX_BAD_SACPZ, IRIS_SACPZDB_FIXES, BAD_SACPZ

%     Version History:
%        Dec.  2, 2009 - initial version
%        May  28, 2010 - handle networks starting with a digit
%        Mar.  5, 2011 - added parallelization
%        Dec. 22, 2011 - convert iris_sacpzdb_howto notes to this script
%        Feb.  3, 2012 - made documentation, removed inert code, use new
%                        form of bad_sacpz_cmplx (avoids crashes)
%        Mar. 15, 2012 - fix parallel verbose
%        Mar. 28, 2013 - turn verbosity back, slight change in output names
%        Mar.  6, 2014 - update for new sacpz struct format
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 11:15 GMT

% todo:

% default path
if(nargin<1 || isempty(path)); path='/opt/seizmo/sacpzdb/sacpz'; end

% build db in matlab
net=xdir(path);
net(strcmp({net.name},'.') | strcmp({net.name},'..'))=[];
verbose=seizmoverbose;
matlabpool(4);
sacpzdb_tmp=cell(numel(net),1);
parfor i=1:numel(net)
    disp(net(i).name);
    seizmoverbose(false);
    sacpzdb_tmp{i}=makesacpzdb([net(i).path net(i).name]);
end
matlabpool close;
seizmoverbose(verbose);

% flatten db & clear paths
disp('FLATTENING DB');
sacpzdb_tmp=sscat(sacpzdb_tmp{:});
sacpzdb_tmp.path=repmat({'./'},numel(sacpzdb_tmp.path),1);

% break into subdbs (for speed)
disp('SPLITTING INTO SUBDBs');
knetwk=sacpzdb_tmp.knetwk;
nets=unique(knetwk);
for i=1:numel(nets)
    % handle networks that begin with a digit
    if(isstrprop(nets{i}(1),'digit'))
        netname=['A_' nets{i}];
    else
        netname=nets{i};
    end
    sacpzdb.(netname)=ssidx(sacpzdb_tmp,strcmpi(knetwk,nets(i)));
end

% save unfixed db (just in case)
disp('FIXING KNOWN BAD RESPONSES & SAVING');
save sacpzdb_unfixed -struct sacpzdb

% apply fixes, mark bad, and save
sacpzdb=fix_bad_sacpz(sacpzdb);
[~,sacpzdb]=bad_sacpz(sacpzdb);
save sacpzdb -struct sacpzdb

end
