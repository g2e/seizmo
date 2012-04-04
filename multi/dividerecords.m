function [data]=dividerecords(varargin)
%DIVIDERECORDS    Divide SEIZMO records
%
%    Usage:    data=dividerecords(data)
%              data=dividerecords(data1,data2)
%              data=dividerecords(data1,data2,...,dataN)
%              data=dividerecords(...,'newhdr',true|false)
%              data=dividerecords(...,'npts',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=dividerecords(...,'ncmp',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=dividerecords(...,'delta','error'|'warn'|'ignore')
%              data=dividerecords(...,'begin','error'|'warn'|'ignore')
%              data=dividerecords(...,'ref','error'|'warn'|'ignore')
%              data=dividerecords(...,'leven','error'|'warn'|'ignore')
%              data=dividerecords(...,'iftype','error'|'warn'|'ignore')
%
%    Description:
%     DATA=DIVIDERECORDS(DATA) will divide records 2+ from record 1 in
%     DATA using the formula:
%         DATA=DATA(1)/DATA(2)/DATA(3)/.../DATA(N)
%     The result is one record with its header fields set to those of the
%     first record.  The header can be set to that of the last record by
%     setting option 'newhdr' to TRUE.  Records should be of the same filetype, be
%     evenly sampled, have the same sample rate, number of points, and
%     timing but these can all be ignored (for better or for worse) by
%     setting options available in BINOPERR to 'ignore'.
%     
%     DATA=DIVIDERECORDS(DATA1,DATA2) will divide the records in DATA1 by
%     DATA2.  If either DATA1 or DATA2 are a single record, the record will
%     be applied to every record of the other dataset.  DATA1 and DATA2
%     must contain the same number of records otherwise.
%     
%     DATA=DIVIDERECORDS(DATA1,DATA2,...,DATAN) divides records in datasets
%     2-N from those in dataset 1.  Every dataset must have the same number
%     of records.  The exception to this rule is single record datasets,
%     which are replicated to match the size of the rest of the datasets.
%     
%     DATA=DIVIDERECORDS(...,'newhdr',true) controls the inheritance of
%     header fields and will set the resultant records' header fields to
%     those of the last record to be divided.  So
%                    DIVIDERECORDS(DATA1,DATA2,'newhdr',true)
%     will produce records with header fields set to those in DATA2.  By
%     default 'newhdr' is set to FALSE which sets the resultant records'
%     header fields to those in DATA1.  If dividing all records in a
%     single dataset, setting 'newhdr' to TRUE will set the resultant
%     record's header equal to the last record's header.  Leaving
%     'newhdr' set to the default FALSE will set the resultant record's
%     header to that of the first record's header.
%     
%     *********************************************************
%     The following options may also be controlled by BINOPERR.
%     *********************************************************
%     
%     DATA=DIVIDERECORDS(...,'npts','error|warn|truncate|pad|ignore') sets
%     the reaction to records with different numbers of points.  If the
%     option is set to 'warn' or 'ignore', the number of points in the
%     records is not altered - which will likely cause an error during the
%     operation.  If the option is set to 'truncate', the number of points
%     in the records will be equal to the record with the least.  Option
%     'pad' will make the records being operated on have number of points
%     equal to that with the most (note that padding is done with zeros).
%     By default 'npts' is set to 'error'.
%     
%     DATA=DIVIDERECORDS(...,'ncmp','error|warn|truncate|pad|ignore') sets
%     the reaction to records with different numbers of components.  If the
%     option is set to 'warn' or 'ignore', the number of components in the
%     records is not altered - which will likely lead to an error.  If the
%     option is set to 'truncate', the number of components in the records
%     being operated on will be equal to that with the least.  Option 'pad'
%     will make the number of components for records in the operation equal
%     to that of the record with the most (note that padding is done with
%     zeros).  By default 'ncmp' is set to 'error'.
%     
%     DATA=DIVIDERECORDS(...,'delta','error|warn|ignore') sets the reaction
%     to records with different sample rates.  If the option is set to
%     'warn' or 'ignore', the records are just operated on point for point
%     (basically ignoring timing).  The resultant records' sample rates are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'delta' is set to 'error'.
%     
%     DATA=DIVIDERECORDS(...,'begin','error|warn|ignore') sets the reaction
%     to records with different begin times.  If the option is set to
%     'warn' or 'ignore', the resultant records' begin times are determined
%     by the parent of their header fields (set by option 'newhdr').  By
%     default 'begin' is set to 'warn'.
%     
%     DATA=DIVIDERECORDS(...,'ref','error|warn|ignore') sets the reaction
%     to records with different reference times.  If the option is set to
%     'warn' or 'ignore', the resultant records' reference times are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'ref' is set to 'warn'.
%     
%     DATA=DIVIDERECORDS(...,'leven','error|warn|ignore') sets the reaction
%     to unevenly sampled records.  If the option is set to 'warn' or
%     'ignore', the records are just operated on point for point (basically
%     ignoring timing).  The resultant records' leven fields are determined
%     by the parent of their header fields (set by option 'newhdr').  By
%     default 'leven' is set to 'error'.
%     
%     DATA=DIVIDERECORDS(...,'iftype','error|warn|ignore') sets the
%     reaction to records of different types.  If the option is set to
%     'warn' or 'ignore', the records are just operated on point for point.  The 
%     resultant records' iftypes are determined by the parent of their
%     header fields (set by option 'newhdr').  By default 'iftype' is set
%     to 'error'.
%     
%    Notes:
%     
%    Header changes: DEPMIN, DEPMAX, DEPMEN,
%     NPTS, NCMP (see option 'npts' and 'ncmp')
%     See option 'newhdr' for inheritance of other header fields.
%
%    Examples:
%     % Plot the time-dependent amplitude ratio between two records:
%     plot1(dividerecords(data(1:2)))
%     
%    See also: MULTIPLYRECORDS, ADDRECORDS, SUBTRACTRECORDS, BINOPERR

%     Version History:
%        June 10, 2008 - Initial Version
%        June 11, 2008 - Full filetype and class support
%        June 20, 2008 - doc update, 'ncmp' option
%        June 29, 2008 - fixed amph2rlim handling, doc update,
%                        .dep and .ind rather than .x and .t
%        Oct.  6, 2008 - doc update, code clean, more checks
%        Nov. 23, 2008 - now just calls RECORDFUN
%        Apr. 23, 2009 - move usage up
%        June 28, 2009 - cleaned up docs for recent changes to RECORDFUN
%        Jan.  6, 2011 - recordfun/multifun rename
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 19:55 GMT

% todo:

% pass on to multifun
data=multifun('/',varargin{:});

end
