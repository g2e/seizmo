function [data,removed]=removedeadrecords(data,option)
%REMOVEDEADRECORDS    Removes constant SEIZMO records
%
%    Usage:    data=removedeadrecords(data)
%              data=removedeadrecords(data,option)
%              [data,removed]=removedeadrecords(...)
%
%    Description:
%     DATA=REMOVEDEADRECORDS(DATA) removes records that have no change in
%     the dependent component.  These can cause problems in analysis and
%     are not worth keeping.  Uses the header fields 'depmin'/'depmax' so
%     that records can be eliminated before reading in the data.
%
%     DATA=REMOVEDEADRECORDS(DATA,OPTION) allows changing how records are
%     determined as dead.  OPTION is a logical that when set true
%     (default), will use the header fields depmin/depmax to look for dead
%     records.  When OPTION is false, the data (in .dep) is used to find
%     dead records.
%
%     [DATA,REMOVED]=REMOVEDEADRECORDS(...) also returns a listing of the
%     indices of the records removed in REMOVED.  These indices are
%     relative to the input dataset, not the output dataset (obviously
%     because the records are no longer in the output dataset!).
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     % Remove dead records before reading in data from current directory:
%     data=readdata(removedeadrecords(readheaders('*')));
%
%    See also: REMOVEMEAN, REMOVETREND, REMOVEPOLYNOMIAL, REMOVESPLINE

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 28, 2008 - SEISCHK support
%        Mar.  4, 2008 - minor doc update
%        Nov. 22, 2008 - update for new name schema (now REMOVEDEADRECORDS)
%        Dec.  8, 2008 - 2nd output lists removed indices
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Jan. 30, 2010 - slimmed the code (no checkheader call)
%        Apr.  1, 2010 - detail message indicates number removed
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix, warn fix
%        Apr.  3, 2012 - minor doc update
%        Jan. 21, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 21, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(1,2,nargin));

% check option
if(nargin==1 || isempty(option))
    option=true;
elseif(~islogical(option) || ~isscalar(option))
    error('seizmo:removedeadrecords:badInput',...
        'OPTION must be true or false!');
end

% detail message
if(seizmoverbose); disp('Removing Dead Record(s)'); end

% proceed by option
if(option)
    % get depmax, depmin
    [depmax,depmin]=getheader(data,'depmax','depmin');
    
    % remove dead records
    removed=((depmax-depmin)==0);
    data(removed)=[];
else
    % check data structure
    error(seizmocheck(data,'dep'));
    
    % get min/max of data
    nrecs=numel(data);
    dmin=nan(nrecs,1); dmax=dmin;
    for i=1:nrecs
        dmin(i)=min(data(i).dep(:));
        dmax(i)=max(data(i).dep(:));
    end
    
    % remove dead records
    removed=((dmax-dmin)==0);
    data(removed)=[];
end

% detail message
if(seizmoverbose)
    disp([' --> Found ' num2str(sum(removed)) ' Dead Record(s)']);
end

end
