function [db]=makesacpzdb(varargin)
%MAKESACPZDB    Makes a SAC PoleZero database using a list of directories
%
%    Usage:    db=makesacpzdb()
%              db=makesacpzdb(dir1,...,dirN)
%
%    Description:
%     DB=MAKESACPZDB() creates a structure from all valid SAC Polezero
%     files in the current directory.  See the Notes sections of
%     PARSE_SACPZ_FILENAME and READSACPZ for file format info.
%
%     SACPZ data structure layout:
%      path   - path to PoleZero file
%      name   - PoleZero file name
%      knetwk - network associated with PoleZero file
%      kstnm  - station associated with PoleZero file
%      kcmpnm - component associated with PoleZero file
%      khole  - stream associated with PoleZero file
%      b      - begin time for PoleZero validity
%      e      - end time for PoleZero validity
%      z      - zeros
%      p      - poles
%      k      - constant
%
%     DB=MAKESACPZDB(DIR1,...,DIRN) builds the structure using directories
%     DIR1 thru DIRN.
%
%    Notes:
%
%    Examples:
%     % To add-in some SAC PoleZero files to the main SACPZ database:
%      % load and flatten the main db
%      sacpzdb_tmp=load('sacpzdb');
%      sacpzdb_tmp=struct2cell(sacpzdb_tmp);
%      sacpzdb_tmp=cat(1,sacpzdb_tmp{:});
%      
%      % create db from personal directory
%      mysacpzdb=makesacpzdb('my/sacpz/dir');
%
%      % concatenate the dbs
%      sacpzdb_tmp=cat(1,sacpzdb_tmp,mysacpzdb)
%
%      % break up db by network
%      knetwk={sacpzdb_tmp.knetwk};
%      nets=unique(knetwk);
%      for i=1:numel(nets)
%          sacpzdb.(nets{i})=sacpzdb_tmp(strcmpi(knetwk,nets(i)));
%      end
%      
%      % save (NOTE: fix the path)
%      save seizmo/response/sacpzdb -struct sacpzdb
%
%    See also: GETSACPZ, READSACPZ, PARSE_SACPZ_FILENAME, WRITESACPZ,
%              APPLYSACPZ, REMOVESACPZ, DB2SACPZ, GENSACPZNAME

%     Version History:
%        Sep. 20, 2009 - initial version
%        Sep. 23, 2009 - added informative output
%        Nov.  2, 2009 - seizmoverbose support, fixed example
%        May  22, 2010 - progress bar only if verbose
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  3, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));
nin=nargin;

% default to current directory if none given
if(nargin==0); varargin{1}='.'; nin=1; end

% check varargin
if(~iscellstr(varargin))
    error('seizmo:makesacpzdb:badInput','DIR must be a string!');
end

% verbosity
verbose=seizmoverbose;

% loop over directories
db=cell(nin,1);
for i=1:nin
    % detail message
    if(verbose)
        disp(['SAC PoleZero Directory: ' varargin{i}])
    end
    
    % get filelist for this directory
    files=xdir(varargin{i});
    filenames={files.name}.';
    
    % detail message
    if(verbose)
        disp(['  Found ' sprintf('%d',numel(files)) ' Files']);
    end
    
    % filter out invalid and parse valid
    [good,knetwk,kstnm,kcmpnm,khole,b,e]=parse_sacpz_filename(filenames);
    files=files(good);
    
    % handle empty case
    if(isempty(files))
        if(verbose)
            disp('  Found 0 Properly Formatted Filenames');
        end
        db{i}([],1)=struct('path',[],'name',[],'knetwk',[],'kstnm',[],...
            'kcmpnm',[],'khole',[],'b',[],'e',[],'z',[],'p',[],'k',[]);
        continue;
    end
    
    % allocate directory sacpzdb
    n=sum(good);
    db{i}=struct('path',{files.path}.','name',{files.name}.',...
        'knetwk',knetwk,'kstnm',kstnm,'kcmpnm',kcmpnm,'khole',khole,...
        'b',mat2cell(b,ones(n,1)),'e',mat2cell(e,ones(n,1)),'z',[],...
        'p',[],'k',[]);
    
    % detail message
    if(verbose)
        disp(['  Found ' sprintf('%d',n) ' Properly Formatted Filenames']);
        print_time_left(0,n);
    end
    
    % now read in the PoleZero info
    good=true(n,1);
    for j=1:n
        try
            [db{i}(j).z,db{i}(j).p,db{i}(j).k]=...
                readsacpz(fullfile(files(j).path,files(j).name));
        catch
            good(j)=false;
        end
        if(verbose); print_time_left(j,n); end
    end
    
    % remove invalid
    db{i}=db{i}(good,1);
    
    % detail message
    if(verbose)
        disp(...
            ['  Read In ' sprintf('%d',sum(good)) ' SAC PoleZero Files!']);
    end
end

% combine directory dbs
db=cat(1,db{:});

end
