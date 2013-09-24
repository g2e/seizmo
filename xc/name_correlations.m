function [xc]=name_correlations(xc)
%NAME_CORRELATIONS    Creates correlogram filenames based on header fields
%
%    Usage:    xc=name_correlations(xc)
%
%    Description:
%     XC=NAME_CORRELATIONS(XC) sets the .name field for all correlograms in
%     SEIZMO dataset XC based on their header fields.  All other records
%     are ignored.
%
%    Notes:
%
%    Examples:
%     % Zero out the master & slave indices and rename:
%     xc=name_correlations(ch(xc,'user0',0,'user1',0));
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS, ISXC,
%              SPLIT_AUTO_CORRELATIONS, HORZ_CORRELATIONS_SETS,
%              NO_REDUNDANT_CORRELATIONS, IS_FULL_MATRIX_OF_CORRELATIONS,
%              LISTFILES, GENNAME, CHANGENAME

%     Version History:
%        Sep.  9, 2013 - initial version
%        Sep. 19, 2013 - properly optimized checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 19, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check structure
error(seizmocheck(xc));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% safely get necessary header info
try
    [m,s,mi,si,mnm,snm]=getheader(xc,...
        'kuser0','kuser1','user0','user1','kt','kname');
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% simple correlogram checking method
% - require KUSER0=MASTER & KUSER1=SLAVE
ok=strcmp(m,'MASTER') & strcmp(s,'SLAVE');

% rename
d=['%0' num2str(fix(log10(max([mi(ok);si(ok)])))+1) 'd'];
name=strcat('CORR_-_MASTER_-_REC',num2str(mi(ok),d),'_-_',mnm(ok,1),'.',...
    mnm(ok,2),'.',mnm(ok,3),'.',mnm(ok,4),'_-_SLAVE_-_REC',...
    num2str(si(ok),d),'_-_',snm(ok,1),'.',snm(ok,2),'.',snm(ok,3),'.',...
    snm(ok,4));
[xc(ok).name]=deal(name{:});

end
