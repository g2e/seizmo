function [db]=gensacpzname(db)
%GENSACPZNAME    Generates SAC Polezero names based on sacpzdb info
%
%    Usage:    sacpzdb=gensacpzname(sacpzdb)
%
%    Description:
%     SACPZDB=GENSACPZNAME(SACPZDB) generates standard SAC PoleZero
%     filenames based on the fields 'knetwk' 'kstnm', 'kcmpnm', 'khole',
%     'b', and 'e'.  These names are written to the 'name' field of the SAC
%     PoleZero db.
%
%    Notes:
%     - Output name has the following form:
%        SAC_PZs_NT_STA_CMP_LL_YYYY.DDD.HH.MM.SS.FFF_YYYY.DDD.HH.MM.SS.FFF
%       where
%        NT   = Network Name (ie IU, XB, etc) - at least 1 char
%        STA  = Station Name (ie T033, CMB, etc) - at least 1 char
%        CMP  = Component Name (ie BHZ, LH1, etc) - at least 1 char
%        LL   = Stream Name (ie __, 01, 02) - 0-2 chars
%        YYYY = Year - 1+ digits
%        DDD  = Julian day - 1+ digits
%        HH   = Hour - 1+ digits
%        MM   = Minute - 1+ digits
%        SS   = Second - 1+ digits
%        FFF  = Fractional seconds - 1+ digits
%
%    Examples:
%     % Make new names and write out SAC PoleZero files:
%     db2sacpz(gensacpzname(mysacpzdb));
%
%    See also: MAKESACPZDB, GETSACPZ, APPLYSACPZ, REMOVESACPZ,
%              READSACPZ_RDSEED, WRITESACPZ_RDSEED, ISSACPZ_RDSEED,
%              READSACPZ, WRITESACPZ, PARSE_SACPZ_FILENAME, FIX_OLD_SACPZ

%     Version History:
%        Sep. 22, 2009 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  3, 2012 - doc update
%        Mar.  6, 2014 - update for new sacpz struct format
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check db
if(~isstruct(db))
    error('seizmo:gensacpzname:badInput','SACPZDB must be a struct!');
end
reqf={'knetwk' 'kstnm' 'kcmpnm' 'khole' 'b' 'e'};
if(any(~ismember(reqf,fieldnames(db))))
    error('seizmo:gensacpzname:badDB',...
        ['DB must have the following fields:\n' ...
        sprintf('%s ',reqf{:})]);
end

if(issacpz_rdseed(db))
    % convert times to strings
    b=[db.b(:,1:4) fix(round(1e4*db.b(:,5))/1e4) ...
        mod(round(1e4*db.b(:,5)),1e4)];
    b=reshape(sprintf('%04d.%03d.%02d.%02d.%02d.%04d',b.'),22,[]).';
    e=[db.e(:,1:4) fix(round(1e5*db.e(:,5))/1e5) ...
        mod(round(1e5*db.e(:,5)),1e5)];
    e=reshape(sprintf('%04d.%03d.%02d.%02d.%02d.%05d',e.'),23,[]).';
    
    % form names
    db.name=strcat({'SAC_PZs_'},db.knetwk,'_',db.kstnm,'_',...
        db.kcmpnm,'_',db.khole,'_',b,'_',e);
else % old db format
    % convert times to strings
    b=cell2mat({db.b}.');
    b=[b(:,1:4) fix(round(1e4*b(:,5))/1e4) mod(round(1e4*b(:,5)),1e4)];
    b=reshape(sprintf('%04d.%03d.%02d.%02d.%02d.%04d',b.'),22,[]).';
    e=cell2mat({db.e}.');
    e=[e(:,1:4) fix(round(1e5*e(:,5))/1e5) mod(round(1e5*e(:,5)),1e5)];
    e=reshape(sprintf('%04d.%03d.%02d.%02d.%02d.%05d',e.'),23,[]).';
    
    % form name
    name=strcat({'SAC_PZs_'},{db.knetwk}.','_',{db.kstnm}.','_',...
        {db.kcmpnm}.','_',{db.khole}.','_',b,'_',e);
    
    % add names into db
    [db.name]=deal(name{:});
end

end
