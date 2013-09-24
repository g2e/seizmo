function [data]=genname(data,style)
%GENNAME    Generic filenames for SEIZMO records
%
%    Usage:    data=genname(data)
%              data=genname(data,style)
%
%    Description:
%     DATA=GENNAME(DATA) generates a generic filename for each
%     record in SEIZMO struct DATA using the KNAME header field group.
%     Care is taken to avoid duplicate filenames in the same dataset but
%     not with existing files or other datasets.  This is useful for simple
%     and short filenames that indicate the station the record corresponds
%     to.  The .name field of each record is replaced by:
%
%         SAC.KNETWK.KSTNM.KHOLE.KCMPNM.IDX
%
%     where IDX is a 2+ digit number starting at 00 and incremented for
%     records with the same KNAME header field set.
%
%     DATA=GENNAME(DATA,STYLE) generates a filename according to a style
%     specified by STYLE.  STYLE may be one of the following:
%      'seizmo' - (default) see above for format
%      'rdseed' - YYYY.DDD.HH.MM.SS.FFFF.KNETWK.KSTNM.KHOLE.KCMPNM.D.SAC
%     The time in the 'rdseed' style is based on the begin time of the
%     record.
%
%    Notes:
%     - The default matches with my perl code naming scheme.
%
%    Examples:
%     Read, set the name, write:
%      writeseizmo(genname(readseizmo('*')));
%
%    See also: CHANGENAME, LISTFILES, FIX_RDSEED_V48, NAME_CORRELATIONS

%     Version History:
%        June 13, 2010 - initial version
%        Oct.  1, 2010 - added style option (with 'rdseed')
%        Sep.  9, 2013 - some see also additions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  9, 2013 at 01:45 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default style
if(nargin==1 || isempty(style)); style='seizmo'; end

% check style
validstyles={'seizmo' 'rdseed'};
if(~isstring(style) || ~any(strcmpi(style,validstyles)))
    error('seizmo:genname:badInput',...
        ['STYLE must be one of the following:\n' ...
        sprintf('''%s'' ',validstyles{:}) '!']);
end

% operate based on style
switch lower(style)
    case 'seizmo'
        % get necessary info
        kname=getheader(data,'kname');

        % combine into a simple string
        kname=upper(kname);
        kname=strcat(kname(:,1),'.',kname(:,2),'.',...
            kname(:,3),'.',kname(:,4));

        % pass through unique for finding repeats
        [uname,idx,idx]=unique(kname);

        % occurance index (starts with zero)
        oi=nan(size(idx));
        for i=1:numel(uname)
            cmps=idx==i;
            oi(cmps)=0:sum(cmps)-1;
        end

        % now finish off the string
        digits=num2str(max(2,log10(max(oi))+1)); % at least 2 digits
        kname=strcat('SAC.',kname,'.',num2str(oi,['%0' digits 'd']));
    case 'rdseed'
        % get necessary info
        [butc,kname]=getheader(data,'b utc','kname');
        butc=cell2mat(butc);
        
        % convert time to string
        datestr=strcat(num2str(butc(:,1),'%d'),...
            '.',num2str(butc(:,2),'%03d'),...
            '.',num2str(butc(:,3),'%02d'),'.',num2str(butc(:,4),'%02d'),...
            '.',num2str(fix(butc(:,5)),'%02d'),...
            '.',num2str(fix(10000*(butc(:,5)-fix(butc(:,5)))),'%04d'));
        
        % fix khole
        kname(strcmp(kname(:,3),'__'),3)={''};
        
        % combine kname into a simple string
        kname=upper(kname);
        kname=strcat(kname(:,1),'.',kname(:,2),'.',...
            kname(:,3),'.',kname(:,4));
        
        % bring it all together
        kname=strcat(datestr,'.',kname,'.D.SAC');
end

% rename
[data.name]=deal(kname{:});

end
