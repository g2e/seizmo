function [ok]=webinstall_extras(dlflag)
%WEBINSTALL_EXTRAS    Install extra SEIZMO components
%
%    Usage:    ok=webinstall_extras
%              ok=webinstall_extras(true)
%
%    Description:
%     OK=WEBINSTALL_EXTRAS downloads & installs some extra files for SEIZMO
%     such as a SAC polezero db, feature data for mapping, & 3D velocity
%     models.  The download is large at 50+ megabytes so make sure you have
%     a good connection or this operation will take a while and you cannot
%     do anything else at the Matlab/Octave prompt while waiting for the
%     files to download.
%
%     OK=WEBINSTALL_EXTRAS(TRUE) forces a redownload of the extra packages.
%
%    Notes:
%
%    Examples:
%     % Reinstall:
%     uninstall_irisws & webinstall_extras(true)
%
%    See also: UNINSTALL_IRISWS, WEBINSTALL_GSHHS, UNINSTALL_GSHHS,
%              WEBINSTALL_MMAP, UNINSTALL_MMAP, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXPORTFIG, WEBINSTALL_NJTBX, UNINSTALL_NJTBX,
%              WEBINSTALL_EXTRAS, UNINSTALL_EXTRAS, WEBINSTALL_TAUP,
%              UNINSTALL_TAUP, INSTALL_SEIZMO, UNINSTALL_SEIZMO

%     Version History:
%        Feb. 21, 2014 - initial version
%        Feb. 27, 2014 - minor doc update
%        Mar.  1, 2014 - iscfmdb not ready yet
%        Mar.  2, 2014 - minor fix to web links
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  2, 2014 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% directory separator
fs=filesep;

% path to seizmo directory
mypath=fileparts(mfilename('fullpath'));

% check path
if(~exist(mypath,'dir'))
    error('seizmo:webinstall_extras:badPath',...
        ['Directory (' mypath ') does not exist!']);
end

% default/check download flag
if(nargin<1 || isempty(dlflag)); dlflag=false; end
if(~islogical(dlflag) || ~isscalar(dlflag))
    error('seizmo:webinstall_extras:badInput',...
        'DLFLAG must be TRUE or FALSE!');
end

% attempt extras install
try
    % go to desired install location
    cwd=pwd;
    cd(mypath);
    
    % current versions
    sacpzdb='seizmo_iris_sacpzdb.zip';     % ~20mb
    sz3dmod='seizmo_3d_models.zip';        % ~20mb
    mapfeat='seizmo_mapping_features.zip'; % ~10mb
    %iscfmdb='seizmo_iscfmdb.zip';          %  ~6mb
    
    % grab files (either locally or remotely)
    url0='http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/';
    %url1='https://s3-us-west-2.amazonaws.com/seizmo/'; % backup currently
    fprintf(' Getting %s\n',sacpzdb);
    if(~dlflag && exist(sacpzdb,'file'))
        if(~exist([mypath fs sacpzdb],'file'))
            copyfile(which(sacpzdb),'.');
        end
    else
        urlwrite([url0 sacpzdb],sacpzdb);
    end
    fprintf(' Getting %s\n',sz3dmod);
    if(~dlflag && exist(sz3dmod,'file'))
        if(~exist([mypath fs sz3dmod],'file'))
            copyfile(which(sz3dmod),'.');
        end
    else
        urlwrite([url0 sz3dmod],sz3dmod);
    end
    fprintf(' Getting %s\n',mapfeat);
    if(~dlflag && exist(mapfeat,'file'))
        if(~exist([mypath fs mapfeat],'file'))
            copyfile(which(mapfeat),'.');
        end
    else
        urlwrite([url0 mapfeat],mapfeat);
    end
    %fprintf(' Getting %s\n',iscfmdb);
    %if(~dlflag && exist(iscfmdb,'file'))
    %    if(~exist([mypath fs iscfmdb],'file'))
    %        copyfile(which(iscfmdb),'.');
    %    end
    %else
    %    urlwrite([url0 iscfmdb],iscfmdb);
    %end
    
    % unpack files
    % - These should land in their appropriate
    %   locations if I packaged them correctly.
    unzip(sacpzdb);
    unzip(sz3dmod);
    unzip(mapfeat);
    %unzip(iscfmdb);
    
    % return
    cd(cwd);
    
    % all good
    ok=true;
catch
    le=lasterror;
    warning(le.identifier,le.message);
    ok=false;
    cd(cwd);
end

end
