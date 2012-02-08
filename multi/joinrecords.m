function [data]=joinrecords(varargin)
%JOINRECORDS    Join SEIZMO records into multiple-component record(s)
%
%    Usage:    data=joinrecords(data)
%              data=joinrecords(data1,data2)
%              data=joinrecords(data1,data2,...,dataN)
%              data=joinrecords(...,'newhdr',true|false)
%              data=joinrecords(...,'npts',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=joinrecords(...,'ncmp',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=joinrecords(...,'delta','error'|'warn'|'ignore')
%              data=joinrecords(...,'begin','error'|'warn'|'ignore')
%              data=joinrecords(...,'ref','error'|'warn'|'ignore')
%              data=joinrecords(...,'leven','error'|'warn'|'ignore')
%              data=joinrecords(...,'iftype','error'|'warn'|'ignore')
%
%    Description:
%     JOINRECORDS(DATA) will horizontally concatenate all records in DATA,
%     returning one record with as many components as all records in DATA
%     combined.  The header fields are set to those of the first record in
%     DATA.  The header can be set to that of the last record by setting
%     option 'newhdr' to TRUE.  Records should be of the same filetype, be
%     evenly sampled, have the same sample rate, number of points, and
%     timing but these can all be ignored (for better or for worse) by
%     setting options available in BINOPERR to 'ignore'.
%     
%     JOINRECORDS(DATA1,DATA2) will join the records in DATA2 to DATA1.  If
%     either DATA1 or DATA2 are a single record the record will be joined
%     with every record of the other dataset.  DATA1 and DATA2 must contain
%     the same number of records otherwise.
%     
%     JOINRECORDS(DATA1,DATA2,...,DATAN) joins records in all N datasets
%     such that the first record of each dataset is joined together and so
%     on.  Therefore every dataset must have the same number of records or
%     be single record datasets, which are replicated to match the size of
%     the rest of the datasets.
%     
%     JOINRECORDS(...,'newhdr',true) controls the inheritance of header
%     fields and will set the resultant records' header fields to those of
%     the last record to be joined on.  So
%                    JOINRECORDS(DATA1,DATA2,'newhdr',true)
%     will produce records with header fields set to those in DATA2.  By
%     default 'newhdr' is set to FALSE which sets the resultant records'
%     header fields to those in DATA1.  If joining all records in a single
%     dataset, setting 'newhdr' to TRUE will set the resultant record's
%     header equal to the last record's header.  Leaving 'newhdr' set to
%     the default FALSE will set the resultant record's header to that of
%     the first record's header.
%     
%     *********************************************************
%     The following options may also be controlled by BINOPERR.
%     *********************************************************
%     
%     JOINRECORDS(...,'npts','error|warn|truncate|pad|ignore') sets the
%     reaction to records with different numbers of points.  If the option
%     is set to 'warn' or 'ignore', the number of points in the records is
%     not altered - which will likely cause an error during the operation.
%     If the option is set to 'truncate', the number of points in the
%     records being operated on will be equal to that with the least.
%     Option 'pad' will make the records being operated on have number of
%     points equal to that with the most (note that padding is done with
%     zeros).  By default 'npts' is set to 'error'.
%     
%     JOINRECORDS(...,'ncmp','error|warn|truncate|pad|ignore') sets the
%     reaction to records with different numbers of components.  If the
%     option is set to 'warn' or 'ignore', the number of components in the
%     records is not altered before the operation.  If the option is set to
%     'truncate', the number of components in the records being operated on
%     will be equal to that with the least.  Option 'pad' will make the
%     number of components for records in the operation equal to that of
%     the record with the most (note that padding is done with zeros).  By
%     default 'ncmp' is set to 'ignore' for JOINRECORDS and should not be
%     changed.
%     
%     JOINRECORDS(...,'delta','error|warn|ignore') sets the reaction to
%     records with different sample rates.  If the option is set to 'warn'
%     or 'ignore', the records are just operated on point for point
%     (basically ignoring timing).  The resultant records' sample rates are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'delta' is set to 'error'.
%     
%     JOINRECORDS(...,'begin','error|warn|ignore') sets the reaction to
%     records with different begin times.  If the option is set to 'warn'
%     or 'ignore', the resultant records' begin times are determined by the
%     parent of their header fields (set by option 'newhdr').  By default
%     'begin' is set to 'warn'.
%     
%     JOINRECORDS(...,'ref','error|warn|ignore') sets the reaction to
%     records with different reference times.  If the option is set to
%     'warn' or 'ignore', the resultant records' reference times are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'ref' is set to 'warn'.
%     
%     JOINRECORDS(...,'leven','error|warn|ignore') sets the reaction to
%     unevenly sampled records.  If the option is set to 'warn' or
%     'ignore', the records are just operated on point for point (basically
%     ignoring timing).  The resultant records' leven fields are determined
%     by the parent of their header fields (set by option 'newhdr').  By
%     default 'leven' is set to 'error'.
%     
%     JOINRECORDS(...,'iftype','error|warn|ignore') sets the reaction to
%     records of different types.  If the option is set to 'warn' or
%     'ignore', the records are just operated on point for point.  The 
%     resultant records' iftypes are determined by the parent of their
%     header fields (set by option 'newhdr').  By default 'iftype' is set
%     to 'error'.
%
%    Notes:
%     - See functions MELD or MULTIFUN to concatenate records time-wise
%       (aka vertically).
%    
%    Header changes: DEPMEN, DEPMIN, DEPMAX,
%     NPTS, NCMP (see options 'npts' and 'ncmp')
%     See option 'newhdr' for inheritance of other header fields.
%
%    Examples:
%     % Separate spectral record components then put them back together:
%     data=joinrecords(keepam(data),keepph(data));
%
%    See also: SPLITRECORDS, MULTIFUN, MELD, ADDRECORDS, SUBTRACTRECORDS,
%              DIVIDERECORDS, MULTIPLYRECORDS, BINOPERR

%     Version History:
%        June 28, 2009 - initial version
%        Jan. 30, 2010 - minor doc update
%        Apr. 25, 2010 - updates ncmp stuff now, doc fix
%        Jan.  6, 2011 - recordfun/multifun rename
%        Feb.  7, 2012 - merge to meld update, doc update
%     
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 19:55 GMT

% todo:

% pass on to multifun
data=multifun(@(x,y)([x y]),'ncmp','ignore',varargin{:});

% fix ncmp
old=checkheader_state(true);
data=checkheader(data,...
    'all','ignore',...
    'inconsistent_dep_ncmp','fix',...
    'invalid_mulcmp_dep','fix');
checkheader_state(old);

end
