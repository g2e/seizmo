function [db]=makesacpzdb(varargin)
%MAKESACPZDB    Makes a SAC PoleZero database using a list of directories
%
%    Usage:    db=makesacpzdb()
%              db=makesacpzdb(dir1,...,dirN)
%
%    Description: DB=MAKESACPZDB() creates a structure from all valid SAC
%     Polezero files in the current directory.  See the Notes sections of
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
%     To add-in some SAC PoleZero files to the main SACPZ database:
%      load sacpzdb
%      mysacpzdb=makesacpzdb('my/sacpz/dir');
%      sacpzdb=cat(1,sacpzdb,mysacpzdb)
%      save seizmo/data/sacpzdb sacpzdb  # NOTE: fix the path
%
%    See also: GETSACPZ, READSACPZ, PARSE_SACPZ_FILENAME, WRITESACPZ,
%              APPLYSACPZ, REMOVESACPZ, DB2SACPZ, GENSACPZNAME

%     Version History:
%        Sep. 20, 2009 - initial version
%        Sep. 23, 2009 - added informative output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 23, 2009 at 21:35 GMT

% todo:

% check nargin
msg=nargchk(0,1,nargin);
if(~isempty(msg)); error(msg); end
nin=nargin;

% default to current directory if none given
if(nargin==0); varargin{1}='.'; nin=1; end

% check varargin
if(~iscellstr(varargin))
    error('seizmo:makesacpzdb:badInput','DIR must be a string!');
end

% loop over directories
db=cell(nin,1);
for i=1:nin
    % status update
    disp(['SAC PoleZero Directory: ' varargin{i}])
    
    % get filelist for this directory
    files=xdir(varargin{i});
    filenames={files.name}.';
    
    % status update
    disp(['  Found ' sprintf('%d',numel(files)) ' Files']);
    
    % filter out invalid and parse valid
    [good,knetwk,kstnm,kcmpnm,khole,b,e]=parse_sacpz_filename(filenames);
    files=files(good);
    
    % handle empty case
    if(isempty(files))
        disp('  Found 0 Properly Formatted Filenames');
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
    
    % status update
    disp(['  Found ' sprintf('%d',n) ' Properly Formatted Filenames']);
    
    % now read in the PoleZero info
    good=true(n,1);
    print_time_left(0,n);
    for j=1:n
        try
            [db{i}(j).z,db{i}(j).p,db{i}(j).k]=...
                readsacpz(fullfile(files(j).path,files(j).name));
        catch
            good(j)=false;
        end
        print_time_left(j,n);
    end
    
    % remove invalid
    db{i}=db{i}(good,1);
    
    % status update
    disp(['  Read In ' sprintf('%d',sum(good)) ' SAC PoleZero Files!'])
end

% combine directory dbs
db=cat(1,db{:});

end
