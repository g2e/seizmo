function []=clean_seizmo(special_string)
%CLEAN_SEIZMO    Deletes Extra SEIZMO Packages
%
%    Usage:    clean_seizmo('YES DELETE!')
%
%    Description:
%     CLEAN_SEIZMO('YES DELETE!') attempts to clean up a SEIZMO install to
%     how it was originally downloaded.  This basically just deletes
%     downloaded/expanded packages.  This could be useful if you are having
%     issues with the external packages.  The string 'YES DELETE!' is
%     required as this function will delete stuff in the seizmo folder.
%
%    Notes:
%
%    Examples:
%     % To do a full reinstall:
%     uninstall_seizmo;
%     clean_seizmo('YES DELETE!');
%     install_seizmo;
%
%    See also: INSTALL_SEIZMO, UNINSTALL_SEIZMO

%     Version History:
%        Mar.  2, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  2, 2014 at 15:25 GMT

% todo:

% root directory
mypath=fileparts(fileparts(mfilename('fullpath')));

% only if as required
if(nargin==1 && isequal(special_string,'YES DELETE!'))
    cwd=pwd;
    cd(mypath);
    delete('*.zip');
    rmdir('export_fig','s');
    rmdir('gshhg','s');
    rmdir('m_map','s');
    rmdir('njtbx','s');
    cd('mattaup/lib');
    delete('TauP-*.jar');
    delete('seisFile-*.jar');
    cd(mypath);
    cd('ws');
    delete('IRIS-WS-*.jar');
    delete('irisFetch.m');
    cd(mypath);
    cd('event');
    delete('globalcmt_*.mat');
    cd(mypath);
    cd('response');
    delete('sacpzdb.mat');
    cd(mypath);
    cd('mapping');
    delete('*.gz');
    delete('*.mat');
    cd(mypath);
    cd('models');
    delete('block.desc');
    delete('CN*.txt');
    delete('crust1.*');
    delete('crust1.*');
    delete('*.mat');
    delete('*.bin');
    cd(mypath);
    delete('CRUST*.mat');
    delete('CUB2.mat');
    delete('DZ04.mat');
    delete('HMSL06*');
    delete('MITP08.mat');
    delete('PRI05.mat');
    delete('S20RTS.mat');
    delete('SAW24B16.mat');
    delete('SB4L18.*');
    delete('TX200*.mat');
    delete('CN*.txt');
    delete('crust1.*');
    delete('block.desc');
    cd(cwd);
else
    error('seizmo:clean_seizmo:badInput',...
        'Are you REALLY sure you want to clean up SEIZMO?');
end

end

