function []=globalcmt_create(overwrite)
%GLOBALCMT_CREATE    Create local copy of GlobalCMT catalog (req internet)
%
%    Usage:    globalcmt_create
%              globalcmt_create(overwrite)
%
%    Description:
%     GLOBALCMT_CREATE will download the GlobalCMT Project's catalog from
%     their website and create a local copy.  This is then accessed by
%     functions like FINDCMT & FINDCMTS.  Use GLOBALCMT_UPDATE as needed to
%     keep the local copy in sync with the website.
%
%     GLOBALCMT_CREATE(OVERWRITE) quietly overwrites an existing local
%     catalog if it exists when OVERWRITE is TRUE.  The default is FALSE.
%
%    Notes:
%     - Needs write permission to SEIZMO directories.
%     - Also updates the cached catalogs under SEIZMO.GLOBALCMT
%
%    Examples:
%     % You will want to update immediately after creating the local copy:
%     globalcmt_create;
%     globalcmt_update;
%
%    See also: GLOBALCMT_UPDATE, READNDK, FINDCMTS, FINDCMT, SETEVENT

%     Version History:
%        Mar.  1, 2012 - initial version
%        Mar. 25, 2013 - no error if cannot write (just warn and move on)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 17:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default overwrite to false
if(nargin<1 || isempty(overwrite)); overwrite=false; end
if(~isscalar(overwrite) || ~islogical(overwrite))
    error('seizmo:globalcmt_create:badInput',...
        'OVERWRITE must be TRUE or FALSE!');
end

% verbosity
verbose=seizmoverbose;

% check for existing copy
catalog='globalcmt_full.mat';
if(exist(catalog,'file'))
    if(~overwrite)
        fprintf('Local CMT Catalog: %s\nFile Exists!\n',catalog);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% SEIZMO global access
global SEIZMO

% location of cmt archives
url='http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/';
files={'jan76_dec10.ndk' ...
    'PRE1976/deep_1962-1976.ndk' 'PRE1976/intdep_1962-1975.ndk'};

% read in archives
cmt=readndk(urlread([url files{1}]),true);
cmt1=readndk(urlread([url files{2}]),true);
cmt2=readndk(urlread([url files{3}]),true);

% remove deep dupes (from 1976)
cmt1=ssidx(cmt1,cmt1.year<1976);

% concatenate archives
cmt=sscat(cmt,cmt1,cmt2);

% remove dupes by name (also sorts by name)
[name,idx]=unique(cmt.name);
if(verbose && numel(name)~=numel(cmt.name))
    dropped=numel(cmt.name)-numel(name);
    disp(['Removed ' num2str(dropped) ' Duplicate CMTs By Name.']);
end
cmt=ssidx(cmt,idx);

% save full
path=fileparts(mfilename('fullpath'));
SEIZMO.GLOBALCMT.FULL=cmt;
try
    if(isoctave)
        save([path filesep 'globalcmt_full.mat'],'-struct','-7','cmt');
    else % matlab
        save([path filesep 'globalcmt_full.mat'],'-struct','cmt');
    end
catch
    le=lasterror;
    warning(le.identifier,le.message);
end

end
