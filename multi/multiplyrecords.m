function [data]=multiplyrecords(varargin)
%MULTIPLYRECORDS    Multiply SEIZMO records
%
%    Usage:    data=multiplyrecords(data)
%              data=multiplyrecords(data1,data2)
%              data=multiplyrecords(data1,data2,...,dataN)
%              data=multiplyrecords(...,'newhdr',true|false)
%              data=multiplyrecords(...,'npts',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=multiplyrecords(...,'ncmp',...
%                  'error'|'warn'|'truncate'|'pad'|'ignore')
%              data=multiplyrecords(...,'delta','error'|'warn'|'ignore')
%              data=multiplyrecords(...,'begin','error'|'warn'|'ignore')
%              data=multiplyrecords(...,'ref','error'|'warn'|'ignore')
%              data=multiplyrecords(...,'leven','error'|'warn'|'ignore')
%              data=multiplyrecords(...,'iftype','error'|'warn'|'ignore')
%
%    Description:
%     DATA=MULTIPLYRECORDS(DATA) will multiply all records in DATA,
%     returning one record with its header fields set to those of the
%     first record.  The header can be set to that of the last record by
%     setting option 'newhdr' to TRUE.  Records should be of the same
%     filetype, be evenly sampled, have the same sample rate, number of
%     points, and timing but these can all be ignored (for better or for
%     worse) by setting options available in BINOPERR to 'ignore'.
%     
%     DATA=MULTIPLYRECORDS(DATA1,DATA2) will multiply the records in DATA1
%     by DATA2.  If either DATA1 or DATA2 are a single record the record
%     will be applied to every record of the other dataset.  DATA1 and
%     DATA2 must contain the same number of records otherwise.
%     
%     DATA=MULTIPLYRECORDS(DATA1,DATA2,...,DATAN) multiplies records in all
%     N datasets.  Every dataset must have the same number of records.  The
%     exception to this rule is single record datasets, which are
%     replicated to match the size of the rest of the datasets.
%     
%     DATA=MULTIPLYRECORDS(...,'newhdr',true) controls the inheritance of
%     header fields and will set the resultant records' header fields to
%     those of the last record to be multiplied.  So
%                    MULTIPLYRECORDS(DATA1,DATA2,'newhdr',true)
%     will produce records with header fields set to those in DATA2.  By
%     default 'newhdr' is set to FALSE which sets the resultant records'
%     header fields to those in DATA1.  If multiplying all records in a
%     single dataset, setting 'newhdr' to TRUE will set the resultant
%     record's header equal to the last record's header.  Leaving
%     'newhdr' set to the default FALSE will set the resultant record's
%     header to that of the first record's header.
%     
%     *********************************************************
%     The following options may also be controlled by BINOPERR.
%     *********************************************************
%     
%     DATA=MULTIPLYRECORDS(...,'npts','error|warn|truncate|pad|ignore')
%     sets the reaction to records with different numbers of points.  If
%     the option is set to 'warn' or 'ignore', the number of points in the
%     records is not altered - which will likely cause an error during the
%     operation.  If the option is set to 'truncate', the number of points
%     in the records will be equal to the record with the least.  Option
%     'pad' will make the records being operated on have number of points
%     equal to that with the most (note that padding is done with zeros).
%     By default 'npts' is set to 'error'.
%     
%     DATA=MULTIPLYRECORDS(...,'ncmp','error|warn|truncate|pad|ignore')
%     sets the reaction to records with different numbers of components.
%     If set to 'warn' or 'ignore', the number of components in the records
%     is not altered - which will likely lead to an error.  If the option
%     is set to 'truncate', the number of components in the records being
%     operated on will be equal to that with the least.  Option 'pad' will
%     make the number of components for records in the operation equal to
%     that of the record with the most (note that padding is done with
%     zeros).  By default 'ncmp' is set to 'error'.
%     
%     DATA=MULTIPLYRECORDS(...,'delta','error|warn|ignore') sets the
%     reaction to records with different sample rates.  If the option is
%     set to 'warn' or 'ignore', the records are just operated on point for
%     point (basically ignoring timing).  The resultant records' sample
%     rates are determined by the parent of their header fields (set by
%     option 'newhdr').  By default 'delta' is set to 'error'.
%     
%     DATA=MULTIPLYRECORDS(...,'begin','error|warn|ignore') sets the
%     reaction to records with different begin times.  If the option is set
%     to 'warn' or 'ignore', the resultant records' begin times are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'begin' is set to 'warn'.
%     
%     DATA=MULTIPLYRECORDS(...,'ref','error|warn|ignore') sets the reaction
%     to records with different reference times.  If the option is set to
%     'warn' or 'ignore', the resultant records' reference times are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'ref' is set to 'warn'.
%     
%     DATA=MULTIPLYRECORDS(...,'leven','error|warn|ignore') sets the
%     reaction to unevenly sampled records.  If the option is set to 'warn'
%     or 'ignore', the records are just operated on point for point
%     (basically ignoring timing).  The resultant records' leven fields are
%     determined by the parent of their header fields (set by option
%     'newhdr').  By default 'leven' is set to 'error'.
%     
%     DATA=MULTIPLYRECORDS(...,'iftype','error|warn|ignore') sets the
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
%     % Convolve a record with itself:
%     idft(multiplyrecords(dft(data([1 1]))))
%     
%    See also: DIVIDERECORDS, ADDRECORDS, SUBTRACTRECORDS, BINOPERR

%     Version History:
%        June 10, 2008 - initial version
%        June 11, 2008 - full filetype and class support
%        June 20, 2008 - doc update, 'ncmp' option
%        Oct.  6, 2008 - doc update, code clean, more checks, added example
%                        fixed amph2rlim handling, .dep and .ind rather 
%                        than .x and .t
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
data=multifun('*',varargin{:});

end
