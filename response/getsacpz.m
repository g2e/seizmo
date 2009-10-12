function [data]=getsacpz(data,varargin)
%GETSACPZ    Finds, reads and adds SAC PoleZero info into a SEIZMO struct
%
%    Usage:    data=getsacpz(data)
%              data=getsacpz(data,sacpzdb)
%              data=getsacpz(data,'dir1',...,'dirN')
%
%    Description: DATA=GETSACPZ(DATA) finds, reads in and adds SAC PoleZero
%     info related to each record into SEIZMO struct DATA.  All info is
%     added under data.misc.sacpz and the format is described in
%     MAKESACPZDB.  This particular case uses the default SAC PoleZero
%     database that comes with SEIZMO.
%
%     DATA=GETSACPZ(DATA,SACPZDB) searches in the alternative database
%     SACPZDB for info relevant to DATA.  See MAKESACPZDB for more info.
%
%     DATA=GETSACPZ(DATA,'DIR1',...,'DIRN') uses directories DIR1 thru DIRN
%     to form a SAC PoleZero database on the fly.  This database is then
%     searched for info relevant to DATA.
%
%    Notes:
%     - warns if no files are found or a record straddles a time boundary
%
%    Examples:
%     Reading in the default full database is slow (>100000 PoleZeros), so
%     to speed operations up in repetitive tasks use the 2nd call type.
%     This requires that you isolate the relevant PoleZero files.  Then
%     just run something similar to the commands below to create a
%     database, save it for later and add the info to the current records:
%      mysacpzdb=makesacpzdb('my/sacpz/dir');
%      save mysacpzdb mysacpzdb
%      data=getsacpz(data,mysacpzdb)
%
%    See also: APPLYSACPZ, REMOVESACPZ, READSACPZ, WRITESACPZ, MAKESACPZDB,
%              PARSE_SACPZ_FILENAME, DB2SACPZ, GENSACPZNAME

%     Version History:
%        Sep. 20, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 20, 2009 at 07:55 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% get sacpzdb
if(nargin==1)
    % load default sacpzdb
    disp('Loading SAC PoleZero Database (May Take Several Minutes)');
    load sacpzdb
elseif(nargin==2 && isstruct(varargin{1}))
    % check struct
    disp('Checking SAC PoleZero Database');
    reqf={'knetwk' 'kstnm' 'kcmpnm' 'khole' 'b' 'e' 'z' 'p' 'k'};
    for i=1:numel(reqf)
        if(~isfield(varargin{1},reqf{i}))
            error('seizmo:getsacpz:badSACPZ',...
                ['SACPZ struct must contain fields:\n' ...
                sprintf('%s ',reqf)]);
        end
    end
    
    % use struct
    sacpzdb=varargin{1};
elseif(iscellstr(varargin))
    % make a sacpzdb using directories given
    disp('Creating SAC PoleZero Database (Gonna Be A While)');
    sacpzdb=makesacpzdb(varargin{:});
else
    error('seizmo:getsacpz:badInputs',...
        'Improper calling of GETSACPZ.  Please try again.');
end

% number of records
nrecs=numel(data);

% get header info
[b,e,knetwk,kstnm,kcmpnm,khole]=getheader(data,'b utc','e utc','knetwk',...
    'kstnm','kcmpnm','khole');

% handle khole goofiness
khole2=khole; badhole=strcmpi(khole,'__');
if(any(badhole)); khole2(badhole)={''}; end

% get db info
dbb=cell2mat({sacpzdb.b}.');
dbe=cell2mat({sacpzdb.e}.');
dbknetwk={sacpzdb.knetwk}.';
dbkstnm={sacpzdb.kstnm}.';
dbkcmpnm={sacpzdb.kcmpnm}.';
dbkhole={sacpzdb.khole}.';

% loop over records
disp('Getting Relevant SAC PoleZeros');
print_time_left(0,nrecs);
for i=1:nrecs
    % set progress bar to overwrite
    redraw=false;
    
    % find sacpz file(s) for this record by name
    % - includes handling khole goofiness
    ok=find(strcmpi(knetwk{i},dbknetwk) & strcmpi(kstnm{i},dbkstnm) ...
        & strcmpi(kcmpnm{i},dbkcmpnm) ...
        & (strcmpi(khole{i},dbkhole) | strcmpi(khole2{i},dbkhole)));
    
    % now find sacpz file(s) for this record by time
    % - find sacpz surrounding record's b & e
    % - warn if overlaps boundary
    okb=timediff(dbb(ok,:),b{i})>0 & timediff(dbe(ok,:),b{i})<0;
    oke=timediff(dbb(ok,:),e{i})>0 & timediff(dbe(ok,:),e{i})<0;
    halfbaked=(okb & ~oke) | (~okb & oke);
    if(any(halfbaked))
        redraw=true;
        warning('seizmo:getsacpz:halfbaked',...
            ['Record: %d\n' ...
            'Record overlaps SAC PoleZero file time boundary!'],i);
    end
    ok=ok(okb & oke);
    
    % warn if no files found
    if(isempty(ok))
        redraw=true;
        warning('seizmo:getsacpz:noGoodSACPZ',...
            'Record: %d\nCould not find a matching SAC PoleZero file!',i);
    else
        % get file with latest b
        [idx,idx]=min(timediff(dbb(ok,:),b{i}));
        
        % assign to data.misc.sacpz
        data(i).misc.sacpz=sacpzdb(ok(idx));
    end
    print_time_left(i,nrecs,redraw);
end

end
