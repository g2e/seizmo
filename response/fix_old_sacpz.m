function [pz]=fix_old_sacpz(pz,butc)
%FIX_OLD_SACPZ    Converts old SEIZMO SAC PoleZero info to new format
%
%    Usage:    pz=fix_old_sacpz(pz)
%              pz=fix_old_sacpz(pz,butc)
%
%    Description:
%     PZ=FIX_OLD_SACPZ(PZ) converts old SEIZMO SAC PoleZero structs to the
%     new RDSEED derived format by putting in place holders for the
%     additional fields.  These fields are not filled with correct info so
%     this is only meant for quick compatibility.
%
%     PZ=FIX_OLD_SACPZ(PZ,BUTC) updates the old SEIZMO SAC PoleZero info
%     using IRIS web services.  BUTC is the start time of the record.  BUTC
%     is necessary as the original polezero time range may not be valid
%     for the corresponding record anymore.  BUTC is expected to be a Nx5
%     array of times as [year jday hour minute seconds].  See the Examples
%     section below to see how to update the polezero info for a SEIZMO
%     dataset using this function.
%
%    Notes:
%
%    Examples:
%     % Update the SAC polezero info for a dataset via IRIS web services:
%     data=getsacpz(data,fix_old_sacpz(getsubfield(data,'misc','sacpz'),...
%         getheader(data,'b utc')));
%
%    See also: ISSACPZ_RDSEED, READSACPZ_RDSEED, WRITESACPZ_RDSEED,
%              READSACPZ, WRITESACPZ, REMOVESACPZ, APPLYSACPZ, MAKESACPZDB,
%              GENSACPZNAME, PARSE_SACPZ_FILENAME, GETSACPZ

%     Version History:
%        Feb. 25, 2014 - initial version
%        Mar.  6, 2014 - update See also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 6, 2014 at 15:05 GMT

% todo:

% number of inputs
error(nargchk(1,2,nargin));

% require old polezero format
if(~issacpz_rdseed(pz,'old'))
    error('seizmo:fix_old_sacpz:badInput',...
        'PZ is not a SAC PoleZero struct in the old format!');
end

% number of polezeros
npz=numel(pz);

% proceed by number of inputs
if(nargin==1)
    % created now
    created=serial2gregorian(now,'doytime');
    
    % add fields
    for i=1:npz
        pz(i).created=created;
        pz(i).description='';
        pz(i).stla=nan;
        pz(i).stlo=nan;
        pz(i).stel=nan;
        pz(i).stdp=nan;
        pz(i).cmpinc=nan;
        pz(i).cmpaz=nan;
        pz(i).sr=nan;
        pz(i).input='';
        pz(i).output='';
        pz(i).insttype='';
        pz(i).instgain=nan;
        pz(i).instgainunits='';
        pz(i).comment='';
        pz(i).sensitivity=nan;
        pz(i).sensitivityunits='';
        pz(i).a0=nan;
    end
    
    % combine
    newpz.path={pz.path}';
    newpz.name={pz.name}';
    newpz.knetwk={pz.knetwk}';
    newpz.kstnm={pz.kstnm}';
    newpz.khole={pz.khole}';
    newpz.kcmpnm={pz.kcmpnm}';
    newpz.created=cell2mat({pz.created}');
    newpz.b=cell2mat({pz.b}');
    newpz.e=cell2mat({pz.e}');
    newpz.description={pz.description}';
    newpz.stla=[pz.stla]';
    newpz.stlo=[pz.stlo]';
    newpz.stel=[pz.stel]';
    newpz.stdp=[pz.stdp]';
    newpz.cmpinc=[pz.cmpinc]';
    newpz.cmpaz=[pz.cmpaz]';
    newpz.sr=[pz.sr]';
    newpz.input={pz.input}';
    newpz.output={pz.output}';
    newpz.insttype={pz.insttype}';
    newpz.instgain=[pz.instgain]';
    newpz.instgainunits={pz.instgainunits}';
    newpz.comment={pz.comment}';
    newpz.sensitivity=[pz.sensitivity]';
    newpz.sensitivityunits={pz.sensitivityunits}';
    newpz.a0=[pz.a0]';
    newpz.z={pz.z}';
    newpz.p={pz.p}';
    newpz.k=[pz.k]';
    pz=newpz;
else
    % check butc
    if(iscell(butc)); butc=cell2mat(butc); end % getheader output style
    if(~isequal(size(butc),[npz 5]) || ~isreal(butc))
        error('seizmo:fix_old_sacpz:badInput',...
            'BUTC must be an Nx5 time array!');
    end
    butc=datestr([doy2cal(butc(:,1:2)) butc(:,3:end)],...
        'yyyy-mm-ddTHH:MM:SS');
    
    % detail message
    verbose=seizmoverbose(false);
    if(verbose)
        disp('Updating Old SAC PoleZero Info');
        print_time_left(0,npz);
    end
    
    % loop over polezeros
    newpz=cell(npz,1);
    for i=1:npz
        % make urls
        url=['http://service.iris.edu/irisws/sacpz/1/query?net=' ...
            pz(i).knetwk '&sta=' pz(i).kstnm '&loc=' pz(i).khole ...
            '&cha=' pz(i).kcmpnm '&time=' butc(i,:)];
    
        % get new polezero info
        try
            newpz{i}=readsacpz_rdseed(urlread(url),true);
        catch exception
            % turn verbose back
            seizmoverbose(verbose);
            
            % rethrow
            throw(exception);
        end
        
        % detail message
        if(verbose); print_time_left(i,npz); end
    end
    
    % turn verbose back
    seizmoverbose(verbose);
    
    % concatenate polezero info
    pz=sscat(newpz{:});
end

end
