function []=db2sacpz(db,varargin)
%DB2SACPZ    Writes out a SAC PoleZero database as PoleZero files
%
%    Usage:    db2sacpz(sacpzdb)
%              db2sacpz(sacpzdb,overwrite)
%
%    Description: DB2SACPZ(SACPZDB) writes out all the entries in the SAC
%     PoleZero database SACPZDB as SAC PoleZero files.  File location and
%     name is determined by the fields 'path' and 'name' in the database.
%     For more info on the fields of a SAC PoleZero database, see the Notes
%     section of makesacpzdb.
%
%     DB2SACPZ(SACPZDB,OVERWRITE) quietly overwrites pre-existing SAC
%     PoleZero files without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.
%
%    Notes:
%     - Attempting to write out the entire sacpzdb included with SEIZMO
%       will take hours to complete and will require ~1GB in disk space!
%     - For your convenience, the 'path' field of sacpzdb is set to '.'
%       (current directory).  This avoids system-dependent path complexity.
%
%    Examples:
%     Write out a network specific (IU) set of SAC PoleZero files:
%      load sacpzdb
%      IUdb=sacpzdb(strcmpi('IU',{sacpzdb.knetwk}));
%      cd my/sacpz/dir
%      db2sacpz(IUdb);
%
%    See also: makesacpzdb, writesacpz, readsacpz, parse_sacpz_filename,
%              applysacpz, removesacpz, getsacpz, gensacpzname

%     Version History:
%        Sep. 22, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 22, 2009 at 06:30 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check db
if(~isstruct(db))
    error('seizmo:db2sacpz:badInput','SACPZDB must be a struct!');
end
reqf={'path' 'name' 'z' 'p' 'k'};
for i=1:numel(reqf)
    if(~isfield(db,reqf{i}))
        error('seizmo:db2sacpz:badDB',...
            ['DB must have the following fields:\n' ...
            sprintf('%s ',reqf{:})]);
    end
end

% loop over entries
for i=1:numel(db)
    % write sacpz
    writesacpz(fullfile(db(i).path,db(i).name),...
        db(i).z,db(i).p,db(i).k,varargin{:});
end

end
