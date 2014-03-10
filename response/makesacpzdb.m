function [db]=makesacpzdb(varargin)
%MAKESACPZDB    Makes a SAC PoleZero database using a list of directories
%
%    Usage:    db=makesacpzdb()
%              db=makesacpzdb(dir1,...,dirN)
%
%    Description:
%     DB=MAKESACPZDB() creates a structure from RDSEED "annotated" SAC
%     PoleZero files in the current directory.  See the Notes sections of
%     READSACPZ for basic file format info and ISSACPZ_RDSEED &
%     READSACPZ_RDSEED for info on the annotation and the returned
%     struct.
%
%     DB=MAKESACPZDB(DIR1,...,DIRN) builds the structure using directories
%     DIR1 thru DIRN.
%
%    Notes:
%     - Non-annotated SAC PoleZero files are no longer allowed by
%       MAKESACPZDB.  To annotate an old-style SAC PoleZero file when you
%       cannot get one from IRIS or where ever else use these functions:
%        READSACPZ + FIX_OLD_SACPZ + WRITESACPZ_RDSEED
%
%    Examples:
%     % To add-in some SAC PoleZero files to the main SACPZ database:
%      % load and flatten the main db
%      sacpzdb_tmp=load('sacpzdb');
%      sacpzdb_tmp=struct2cell(sacpzdb_tmp);
%      sacpzdb_tmp=sscat(sacpzdb_tmp{:});
%      
%      % create db from personal directory
%      mysacpzdb=makesacpzdb('my/sacpz/dir');
%
%      % concatenate the dbs
%      sacpzdb_tmp=sscat(sacpzdb_tmp,mysacpzdb)
%
%      % break up db by network
%      knetwk=sacpzdb_tmp.knetwk;
%      nets=unique(knetwk);
%      for i=1:numel(nets)
%          sacpzdb.(nets{i})=sacpzdb_tmp(strcmpi(knetwk,nets(i)));
%      end
%      
%      % save (NOTE: fix the path)
%      save seizmo/response/sacpzdb -struct sacpzdb
%
%    See also: GETSACPZ, READSACPZ, PARSE_SACPZ_FILENAME, WRITESACPZ,
%              APPLYSACPZ, REMOVESACPZ, GENSACPZNAME, FIX_OLD_SACPZ,
%              READSACPZ_RDSEED, WRITESACPZ_RDSEED, ISSACPZ_RDSEED

%     Version History:
%        Sep. 20, 2009 - initial version
%        Sep. 23, 2009 - added informative output
%        Nov.  2, 2009 - seizmoverbose support, fixed example
%        May  22, 2010 - progress bar only if verbose
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  3, 2012 - doc update
%        Mar.  6, 2014 - update for new sacpz struct format
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,inf,nargin));
nin=nargin;

% default to current directory if none given
if(nargin==0); varargin{1}='.'; nin=1; end

% check varargin
if(~iscellstr(varargin))
    error('seizmo:makesacpzdb:badInput','DIR must be a string!');
end

% verbosity
verbose=seizmoverbose(false);

% loop over directories
db=cell(nin,1);
good=true(nin,1);
for i=1:nin
    % detail message
    if(verbose); disp(['SAC PoleZero Directory: ' varargin{i}]); end
    
    % get filelist for this directory
    files=xdir(varargin{i});
    files(strcmp({files.name},'.') | strcmp({files.name},'..'))=[];
    files=strcat({files.path}',{files.name}');
    
    % detail message
    if(verbose)
        disp(['  Found ' sprintf('%d',numel(files)) ' File(s)']);
    end
    
    % % now read in the PoleZero info
    try
        db{i}=readsacpz_rdseed(files{:});
    catch
        good(i)=false;
        le=lasterror;
        warning(le.identifier,le.message);
        continue;
    end
    
    % detail message
    if(verbose)
        disp(['  Read In ' sprintf('%d',numel(db{i}.k)) ...
            ' SAC PoleZeros!']);
    end
end

% combine directory dbs
if(any(good))
    db=sscat(db{good});
else
    error('seizmo:makesacpzdb:badPZFiles',...
        'No Good PoleZero Files Found!');
end

% restore verbosity
seizmoverbose(verbose);

end
