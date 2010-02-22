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
%     uses the default SAC PoleZero database that comes with SEIZMO.  See
%     the Notes section below for header field requirements!
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
%       for each record, the following fields must be set appropriately:
%        KNETWK, KSTNM, KHOLE, KCMPNM
%        NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC, B, E
%     - you may mix the last three usage forms (note that structs are put
%       in front, followed by sacpzs from dirs, and finally dbs on disk)
%     - warns if no files are found or a record straddles a time boundary
%     - If you have entered the necessary header fields and GETSACPZ does
%       not find a match there are several possibilities:
%        1) the instrument response is not at IRIS
%        2) the instrument response is not in a Poles & Zeros format
%           (this is true for ~25% of the responses at IRIS)
%        3) I haven't updated the sacpzdb recently enough to include your
%           site's new response info.  I update after about 3 months.
%
%    Header changes: see CHECKHEADER
%
%    Examples:
%     Working with the full IRIS database (default) is somewhat slow
%     (with >150000 PoleZeros what do you expect?).  If you can isolate the
%     relevant PoleZero files for your dataset(s) then this operation will
%     be significantly faster.  Run something similar to the commands below
%     to create your custom database, save it for later and add the info to
%     a dataset:
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
%        Nov.  2, 2009 - seizmoverbose support, network specific searching,
%                        fixed example to reflect improvements
%        Jan. 26, 2010 - do not error if no matching polezero found in db
%        Jan. 30, 2010 - fixes for checking state function
%        Feb.  3, 2010 - versioninfo caching
%        Feb. 16, 2010 - require header fields are set
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2010 at 02:15 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure & header
data=checkheader(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt sac polezero lookup
try
    % get undef values
    [h,vi]=versioninfo(data);
    sundef=[h(vi).undef];
    nundef=[sundef.ntype]';
    sundef={sundef.stype}';
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % get header info
    [b,e,knetwk,kstnm,kcmpnm,khole]=getheader(data,'b utc','e utc',...
        'knetwk','kstnm','kcmpnm','khole');
    
    % uncell b/e
    b=cell2mat(b);
    e=cell2mat(e);
    
    % require all fields to be defined
    % - all name fields should not be undefined
    badname=strcmpi(sundef,knetwk) | strcmpi(sundef,kstnm) ...
        | strcmpi(sundef,kcmpnm) | strcmpi(sundef,khole);
    if(any(badname))
        error('seizmo:getsacpz:badName',...
            ['Record(s):\n' sprintf('%d ',find(badname)) ...
            '\nKNETWK, KSTNM, KHOLE, and/or KCMPNM fields not set!']);
    end
    % - all time fields should not be undef, nan, inf
    badtime=sum(nundef(:,ones(5,1))==b | isnan(b) | isinf(b) ...
        | nundef(:,ones(5,1))==e | isnan(e) | isinf(e),2)>0;
    if(any(badtime))
        error('seizmo:getsacpz:badName',...
            ['Record(s):\n' sprintf('%d ',find(badtime)) ...
            '\nB, E, and/or NZ* header fields not set!']);
    end
    
    % handle khole goofiness
    khole2=khole; badhole=strcmpi(khole,'__');
    if(any(badhole)); khole2(badhole)={''}; end
    
    % get sacpzdb
    if(nargin==1)
        % networks of dataset
        NET=upper(knetwk);
        nets=unique(NET);
        
        % detail message
        if(verbose)
            disp('Loading SAC PoleZero Databases For Networks:');
            disp(['  ' sprintf('%s ',nets{:})]);
            disp('This May Take A Minute...');
        end
        
        % load IRIS SAC PoleZero database
        sacpzdb=load('sacpzdb',nets{:});
        if(~isequal(fieldnames(sacpzdb),nets))
            badnets=nets(~ismember(nets,fieldnames(sacpzdb)));
            warning('seizmo:getsacpz:badKNETWK',...
                ['These networks have no response info ' ...
                'in the IRIS database:\n' sprintf('%s ',badnets{:})]);
            nets=nets(ismember(nets,fieldnames(sacpzdb)));
        end
        
        % get db info
        for i=1:numel(nets)
            db.(nets{i}).b=cell2mat({sacpzdb.(nets{i}).b}.');
            db.(nets{i}).e=cell2mat({sacpzdb.(nets{i}).e}.');
            db.(nets{i}).knetwk={sacpzdb.(nets{i}).knetwk}.';
            db.(nets{i}).kstnm={sacpzdb.(nets{i}).kstnm}.';
            db.(nets{i}).kcmpnm={sacpzdb.(nets{i}).kcmpnm}.';
            db.(nets{i}).khole={sacpzdb.(nets{i}).khole}.';
        end
    else
        % detail message
        if(verbose)
            disp('Creating Custom SAC PoleZero Database');
        end
        
        % get structs
        valid=false;
        structs=cellfun('isclass',varargin,'struct');
        try
            sacpzdb.CUSTOM=cat(1,varargin{structs});
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
                sacpzdb.CUSTOM=[sacpzdb.CUSTOM; ...
                    makesacpzdb(varargin{dirs})];
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
                    sacpzdb.CUSTOM=cat(1,sacpzdb.CUSTOM,tmp{:});
                catch
                    error('seizmo:getsacpz:badInputs',...
                        'Alternate SACPZdbs are not concatenateable!');
                end
            end
        end
        
        % check sacpzdb
        if(~valid)
            if(verbose)
                disp('Checking Custom SAC PoleZero Database');
            end
            reqf={'knetwk' 'kstnm' 'kcmpnm' 'khole' 'b' 'e' 'z' 'p' 'k'};
            for i=1:numel(reqf)
                if(~isfield(sacpzdb.CUSTOM,reqf{i}))
                    error('seizmo:getsacpz:badSACPZ',...
                        ['SACPZDBs must contain fields:\n' ...
                        sprintf('%s ',reqf)]);
                end
            end
        end
        
        % get db info
        db.CUSTOM.b=cell2mat({sacpzdb.CUSTOM.b}.');
        db.CUSTOM.e=cell2mat({sacpzdb.CUSTOM.e}.');
        db.CUSTOM.knetwk={sacpzdb.CUSTOM.knetwk}.';
        db.CUSTOM.kstnm={sacpzdb.CUSTOM.kstnm}.';
        db.CUSTOM.kcmpnm={sacpzdb.CUSTOM.kcmpnm}.';
        db.CUSTOM.khole={sacpzdb.CUSTOM.khole}.';
        NET(1:nrecs,1)={'CUSTOM'};
        nets={'CUSTOM'};
    end
    
    % detail message
    if(verbose)
        disp('Getting Relevant SAC PoleZero(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    for i=1:nrecs
        % set progress bar to overwrite
        redraw=false;
        
        % skip records with network not in db
        if(ismember(NET{i},nets))
            % find sacpz file(s) for this record by name
            % - includes handling khole goofiness
            ok=find(strcmpi(knetwk{i},db.(NET{i}).knetwk) ...
                & strcmpi(kstnm{i},db.(NET{i}).kstnm) ...
                & strcmpi(kcmpnm{i},db.(NET{i}).kcmpnm) ...
                & (strcmpi(khole{i},db.(NET{i}).khole) ...
                | strcmpi(khole2{i},db.(NET{i}).khole)));

            % now find sacpz file(s) for this record by time
            % - find sacpz surrounding record's b & e
            % - warn if overlaps boundary
            okb=timediff(db.(NET{i}).b(ok,:),b(i,:))>0 ...
                & timediff(db.(NET{i}).e(ok,:),b(i,:))<0;
            oke=timediff(db.(NET{i}).b(ok,:),e(i,:))>0 ...
                & timediff(db.(NET{i}).e(ok,:),e(i,:))<0;
            halfbaked=(okb & ~oke) | (~okb & oke);
            if(any(halfbaked))
                redraw=true;
                warning('seizmo:getsacpz:halfbaked',...
                    ['Record: %d\n' ...
                    'Record overlaps SAC PoleZero file time boundary!'],i);
            end
            ok=ok(okb & oke);
        else
            % no polezero in db
            ok=[];
        end
        
        % warn if no files found
        if(isempty(ok))
            data(i).misc.has_sacpz=false;
            data(i).misc.sacpz=[];
            redraw=true;
            warning('seizmo:getsacpz:noGoodSACPZ',['Record: %d\n' ...
                'Could not find a matching SAC PoleZero file!'],i);
        else
            % get file with latest b
            [idx,idx]=min(timediff(db.(NET{i}).b(ok,:),b(i,:)));
            
            % assign to data.misc.sacpz
            data(i).misc.has_sacpz=true;
            data(i).misc.sacpz=sacpzdb.(NET{i})(ok(idx));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs,redraw); end
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
