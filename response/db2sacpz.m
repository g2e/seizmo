function []=db2sacpz(db,varargin)
%DB2SACPZ    Writes out a SAC PoleZero database as PoleZero files
%
%    Usage:    db2sacpz(sacpzdb)
%              db2sacpz(sacpzdb,overwrite)
%
%    Description:
%     DB2SACPZ(SACPZDB) writes out all the entries in the SAC PoleZero
%     database SACPZDB as SAC PoleZero files.  File location and name is
%     determined by the fields 'path' and 'name' in the database.  For more
%     info on the fields of a SAC PoleZero database, see the Notes section
%     of MAKESACPZDB.
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
%     % Write out a network specific (IU) set of SAC PoleZero files:
%     db=load('sacpzdb','IU');
%     cd my/sacpz/dir
%     db2sacpz(db.IU);
%
%    See also: MAKESACPZDB, WRITESACPZ, READSACPZ, PARSE_SACPZ_FILENAME,
%              APPLYSACPZ, REMOVESACPZ, GETSACPZ, GENSACPZNAME

%     Version History:
%        Sep. 22, 2009 - initial version
%        Feb.  3, 2010 - updated example
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  3, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

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

% verbosity
verbose=seizmoverbose;

% number of pz
npz=numel(db);

% detail message
if(verbose)
    disp('Writing SAC PoleZero Database as Individual File(s)');
    print_time_left(0,npz);
end

% loop over entries
for i=1:npz
    % write sacpz
    writesacpz(fullfile(db(i).path,db(i).name),...
        db(i).z,db(i).p,db(i).k,varargin{:});
    % detail message
    if(verbose); print_time_left(i,npz); end
end

end
