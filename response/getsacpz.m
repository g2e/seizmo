function [data]=getsacpz(data,varargin)
%GETSACPZ    Finds, reads and adds SAC PoleZero info into a SEIZMO struct
%
%    Usage:    data=getsacpz(data)
%              data=getsacpz(data,sacpzdb1,...,sacpzdbN)
%              data=getsacpz(data,'db1',...,'dbN')
%              data=getsacpz(data,'dir1',...,'dirN')
%
%    Description: DATA=GETSACPZ(DATA) finds, reads in and adds SAC PoleZero
%     info related to each record in SEIZMO struct DATA.  All info is added
%     under DATA.misc.sacpz and the format is described in MAKESACPZDB.
%     Records without SAC PoleZero info have the field DATA.misc.has_sacpz
%     set to FALSE (otherwise it is set to TRUE).  This particular case
%     uses the default SAC PoleZero database that comes with SEIZMO.
%
%     DATA=GETSACPZ(DATA,SACPZDB1,...,SACPZDBN) searches in the alternative
%     databases SACPZDB1 thru SACPZDBN for info relevant to DATA.  See
%     MAKESACPZDB for more info.
%
%     DATA=GETSACPZ(DATA,'DB1',...,'DBN') loads and then searches the
%     alternative databases given by locations DB1 thru DB2.
%
%     DATA=GETSACPZ(DATA,'DIR1',...,'DIRN') uses directories DIR1 thru DIRN
%     to form a SAC PoleZero database on the fly.  This database is then
%     searched for info relevant to DATA.
%
%    Notes:
%     - In order for GETSACPZ to identify the appropriate SAC PoleZero file
%       for each record, the following fields must be set correctly:
%        KNETWK, KSTNM, KHOLE, KCMPNM
%        NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC, B, E
%     - you may mix the last three usage forms (note that structs are put
%       in front, followed by sacpzs from dirs, and finally dbs on disk)
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
%        Oct. 19, 2009 - allow multiple structs, path/to/db and sac/pz/dir,
%                        checks data and properly handles errors, default
%                        sacpzdb now allows individual network access
%        Oct. 22, 2009 - add info about required header fields
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 22, 2009 at 00:20 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
    % number of records
    nrecs=numel(data);
    
    % get header info
    [b,e,knetwk,kstnm,kcmpnm,khole]=getheader(data,'b utc','e utc',...
        'knetwk','kstnm','kcmpnm','khole');
    
    % handle khole goofiness
    khole2=khole; badhole=strcmpi(khole,'__');
    if(any(badhole)); khole2(badhole)={''}; end
    
    % get sacpzdb
    if(nargin==1)
        % load default sacpzdb
        tmp=unique(knetwk);
        disp('Loading SAC PoleZero Databases For Networks:');
        disp(['  ' sprintf('%s ',tmp{:})]);
        disp('  This May Take A Minute');
        tmp=load('sacpzdb',tmp{:});
        tmp=struct2cell(tmp);
        sacpzdb=cat(1,tmp{:});
    else
        disp('Creating SAC PoleZero Database (May Take A While)');
        % get structs
        valid=false;
        structs=cellfun('isclass',varargin,'struct');
        try
            sacpzdb=cat(1,varargin{structs});
        catch
            error('seizmo:getsacpz:badInputs',...
                'Alternate SACPZdbs are not concatenateable!');
        end
        varargin(structs)=[];

        % rest should be strings
        if(~iscellstr(varargin))
            error('seizmo:getsacpz:badInputs',...
                'Improper calling of GETSACPZ.  Please try again.');
        end

        % get dirs
        ndirs=numel(varargin);
        dirs=false(1,ndirs);
        for i=1:ndirs
            dirs(i)=isdir(varargin{i});
        end
        if(any(dirs))
            try
                sacpzdb=[sacpzdb; makesacpzdb(varargin{dirs})];
                valid=true;
            catch
                error('seizmo:getsacpz:badInputs',...
                    'Alternate SACPZdbs are not concatenateable!');
            end
        end
        varargin(dirs)=[];

        % load rest
        ndb=numel(varargin);
        if(ndb)
            for i=1:ndb
                tmp=load(varargin{i});
                tmp=struct2cell(tmp);
                try
                    sacpzdb=cat(1,sacpzdb,tmp{:});
                catch
                    error('seizmo:getsacpz:badInputs',...
                        'Alternate SACPZdbs are not concatenateable!');
                end
            end
        end

        % check sacpzdb
        if(~valid)
            disp('Checking SAC PoleZero Database');
            reqf={'knetwk' 'kstnm' 'kcmpnm' 'khole' 'b' 'e' 'z' 'p' 'k'};
            for i=1:numel(reqf)
                if(~isfield(varargin{1},reqf{i}))
                    error('seizmo:getsacpz:badSACPZ',...
                        ['SACPZ struct must contain fields:\n' ...
                        sprintf('%s ',reqf)]);
                end
            end
        end
    end
    
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
            data(i).misc.has_sacpz=false;
            redraw=true;
            warning('seizmo:getsacpz:noGoodSACPZ',['Record: %d\n' ...
                'Could not find a matching SAC PoleZero file!'],i);
        else
            % get file with latest b
            [idx,idx]=min(timediff(dbb(ok,:),b{i}));
            
            % assign to data.misc.sacpz
            data(i).misc.has_sacpz=true;
            data(i).misc.sacpz=sacpzdb(ok(idx));
        end
        print_time_left(i,nrecs,redraw);
    end
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
