function [data,pz]=removesacpz(data,varargin)
%REMOVESACPZ    Removes SAC PoleZero responses from SEIZMO records
%
%    Usage:    [dataout,good]=removesacpz(datain)
%              [...]=removesacpz(...,'freqlimits',[f1 f2 f3 f4],...)
%              [...]=removesacpz(...,'units',UNITS,...)
%              [...]=removesacpz(...,'idep',IDEP,...)
%              [...]=removesacpz(...,'h2o',H2O,...)
%
%    Description: [DATAOUT,GOOD]=REMOVESACPZ(DATAIN) removes the instrument
%     response from records in DATAIN based on the associated SAC PoleZero
%     info.  If this info has not been added using GETSACPZ, then GETSACPZ
%     is called.  If any records in DATA do not have any associated SAC
%     PoleZero info (as is set in the .misc.has_sacpz struct field by
%     GETSACPZ), the records are not returned in DATAOUT.  The secondary
%     output, GOOD, is a logical array indicating the records in DATAIN
%     that had SAC PoleZero info (.misc.has_sacpz set TRUE).
%
%     [...]=REMOVESACPZ(...,'FREQLIMITS',[F1 F2 F3 F4],...) applies a
%     lowpass and a highpass taper that limits the spectrum of the
%     deconvolved records.  F1 and F2 give the highpass taper limits while
%     F3 and F4 specify the lowpass taper limits.  The highpass taper is
%     zero below F1 and unity above F2.  The lowpass taper is zero above F4
%     and unity below F3.  The tapers are cosine tapers applied in the
%     spectral domain.  This is an acausal filter and should not be used if
%     you want to preserve seismic phase onsets.  Use option H2O to
%     stabilize the deconvolution without the frequency domain tapering.
%     Default FREQLIMITS is [-1 -1 2*NYQ 2*NYQ] where NYQ is the nyquist
%     frequency for any particularly record.  The defaults will not apply
%     any tapering to the records.  Note that FREQLIMITS may be specified
%     as F1 or [F1 F2] or [F1 F2 F3].  In these cases the remaining unset
%     values are set to the defaults.
%
%     [...]=REMOVESACPZ(...,'UNITS',UNITS,...) sets the output units.  It
%     is expected that SAC Polezero files will give displacement in meters
%     when the response is removed (note that a meters to nanometers
%     conversion is done before records are output).  Setting UNITS to
%     'velocity' will convert the records to velocity based upon this
%     assumption.  Valid strings are:
%      'none', 'dis', 'vel', and 'acc'.
%     Using any of the last three will automatically update the IDEP values
%     appropriately (next option).  Default UNITS is 'none'.
%
%     [...]=REMOVESACPZ(...,'IDEP',IDEP,...) sets the output dependent
%     component label stored in the header field 'idep' to IDEP.  Typically
%     this is 'idisp', 'ivel', 'iacc', or 'iunkn'.  May be a cell array of
%     strings to set each record separately.  The default IDEP is 'idisp'.
%
%     [...]=REMOVESACPZ(...,'H2O',H2O,...) sets the waterlevel factor for
%     the spectral division.  H2O by default is 0, which has no effect.
%     Adjusting this value (it must be a positive real value) does
%     stabilize the deconvolution and is particularly useful for cases
%     where using the FREQLIMITS tapering is not an option.  Typical values
%     are in the range 0.001 to 0.1.  H2O may be an array of values to set
%     the waterlevel for each record separately.
%
%    Notes:
%     - SAC PoleZero info should be set to convert machine units to
%       displacement in meters
%     - Output by default is displacement in nanometers (1e-9 meters)
%     - the SCALE field is set to 1 (matches SAC)
%     - In order for GETSACPZ to identify the appropriate SAC PoleZero file
%       for each record, the following fields must be set correctly:
%        KNETWK, KSTNM, KHOLE, KCMPNM
%        NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC, B, E
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, IDEP, SCALE
%
%    Examples:
%     Remove the instrument response for a dataset, converting to velocity
%     and ignoring all periods at >200s while tapering periods from 200 to
%     150 seconds:
%      data=removesacpz(data,'units','vel','freqlimits',[1/200 1/150]);
%
%    See also: APPLYSACPZ, GETSACPZ, DECONVOLVE

%     Version History:
%        Oct. 22, 2009 - initial version
%        Oct. 30, 2009 - added informative output on error
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 30, 2009 at 18:50 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% import SEIZMO info
global SEIZMO

% attempt rest
try
    % attempt to access has_sacpz
    try
        pz=getsubfield([data.misc],'has_sacpz');
    catch
        data=getsacpz(data);
        pz=getsubfield([data.misc],'has_sacpz');
    end
    
    % only use those with polezero info
    data=data(pz);
    
    % number of records
    nrecs=numel(data);
    
    % get header info
    leven=getlgc(data,'leven');
    iftype=getenumid(data,'iftype');
    [npts,ncmp,delta,e]=getheader(data,'npts','ncmp','delta','e');
    
    % get spectral
    rlim=strcmpi(iftype,'irlim');
    amph=strcmpi(iftype,'iamph');
    
    % get nyquist frequency
    nyq=1./(2*delta);
    if(any(rlim | amph)); nyq(rlim | amph)=e(rlim | amph); end
    
    % cannot do xyz records
    if(any(strcmpi(iftype,'ixyz')))
        error('seizmo:removesacpz:badIFTYPE',...
            'Illegal operation on XYZ record(s)!');
    end
    
    % cannot do unevenly sampled records
    if(any(strcmpi(leven,'false')))
        error('seizmo:removesacpz:badLEVEN',...
            'Invalid operation on unevenly sampled records!');
    end
    
    % valid values for strings
    valid.UNITS={'n' 'none' 'd' 'dis' 'disp' 'displacement' ...
        'v' 'vel' 'velo' 'velocity' 'a' 'acc' 'accel' 'acceleration'};
    
    % default options
    option.FREQLIMITS=... % set to be outside range that matters
        [-1*ones(nrecs,2) 2*nyq(:,[1 1])]; % frequency limits, Hz
    option.UNITS='none'; % none/dis/vel/acc
    option.IDEP='idisp'; % iunkn/idisp/ivel/iacc/ivolt
    option.H2O=zeros(nrecs,1); % 0+
    flim=option.FREQLIMITS; % backup (for partial freqlimits entries)
    
    % get options from SEIZMO global
    ME=upper(mfilename);
    try
        fields=fieldnames(SEIZMO.(ME));
        for i=1:numel(fields)
            if(~isempty(SEIZMO.(ME).(fields{i})))
                option.(fields{i})=SEIZMO.(ME).(fields{i});
            end
        end
    catch
    end
    
    % get options from command line
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:removesacpz:badInput',...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % check options
    fields=fieldnames(option);
    for i=1:numel(fields)
        value=option.(fields{i});
        
        % which option
        switch lower(fields{i})
            case 'freqlimits'
                % assure real and correct size
                if(isreal(value) && any(size(value,1)==[1 nrecs]) ...
                        && size(value,2)<=4 && ndims(value)==2)
                    if(size(value,1)==1)
                        option.FREQLIMITS=[value(ones(nrecs,1),:) ...
                            flim(:,(size(value,2)+1):4)];
                    else
                        option.FREQLIMITS=[value ...
                            flim(:,(size(value,2)+1):4)];
                    end
                else
                    error('seizmo:removesacpz:badInput',...
                        'FREQLIMITS must be [F1 F2 F3 F4]!');
                end
            case 'units'
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) ...
                        || ~any(numel(value)==[1 nrecs]) ...
                        || any(~ismember(value,valid.UNITS)))
                    error('seizmo:removesacpz:badInput',...
                        'UNITS has wrong class, size, and/or value(s)!');
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                option.UNITS=value;
                
                % now adjust idep
                d=ismember(value,{'d' 'dis' 'disp' 'displacement'});
                v=ismember(value,{'v' 'vel' 'velo' 'velocity'});
                a=ismember(value,{'a' 'acc' 'accel' 'acceleration'});
                value=option.IDEP;
                if(ischar(value)); value=cellstr(value); end
                if(iscellstr(value))
                    if(isscalar(value)); value=value(ones(nrecs,1),1); end
                    value(d)={'idisp'};
                    value(v)={'ivel'};
                    value(a)={'iacc'};
                    option.IDEP=value;
                end
            case 'idep'
                if(ischar(value)); value=cellstr(value); end
                if(~any(numel(value)==[1 nrecs]))
                    error('seizmo:removesacpz:badInput',...
                        ['IDEP must be scalar or have\n' ...
                        '1 element per record in DATA!']);
                end
                option.IDEP=value;
            case 'h2o'
                if(~isreal(value) ...
                        || (~isscalar(value) && numel(value)~=nrecs) ...
                        || any(value<0))
                    error('seizmo:removesacpz:badInput',...
                        'H2O must be a real positive scalar or array!');
                end
            otherwise
                error('seizmo:removesacpz:badInput',...
                    'Unknown option: %s !',varargin{i});
        end
    end
    
    % sort freqlimits
    flim=sort(option.FREQLIMITS,2);
    
    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep)); continue; end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % convert to complex spectra
        if(amph(i))
            nspts=npts(i);
            sdelta=delta(i);
            tmp=data(i).dep(:,1:2:end).*exp(j*data(i).dep(:,2:2:end));
        elseif(rlim(i))
            nspts=npts(i);
            sdelta=delta(i);
            tmp=complex(data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
        else
            nspts=2^(nextpow2(npts(i))+1);
            sdelta=2*nyq(i)./nspts;
            tmp=fft(data(i).dep,nspts,1);
        end
        
        % get limited frequency range
        freq=abs([linspace(0,nyq(i),nspts/2+1) ...
            linspace(-nyq(i)+sdelta,-sdelta,nspts/2-1)]);
        good=freq>=flim(i,1) & freq<=flim(i,4);
        
        % taper
        taper1=taperfun('hann',freq,flim(i,1:2)).';
        taper2=taperfun('hann',nyq(i)-freq,nyq(i)-flim(i,3:4)).';
        tmp=tmp.*taper1(:,ones(ncmp(1),1)).*taper2(:,ones(ncmp(1),1));
        
        % convert zpk to fap
        switch lower(option.UNITS{i})
            case {'n' 'none' 'd' 'dis' 'disp' 'displacement'}
                wpow=0;
            case {'v' 'vel' 'velo' 'velocity'}
                wpow=1;
            case {'a' 'acc' 'accel' 'acceleration'}
                wpow=2;
        end
        [a,p]=zpk2ap(freq,data(i).misc.sacpz.z,data(i).misc.sacpz.p,...
            data(i).misc.sacpz.k,wpow);
        h=((a+option.H2O(i)).*exp(j*p)).';
        
        % remove response (over limited freqrange)
        % - multiply by 1e9 to account for SAC PoleZero in meters
        tmp(good,:)=1e9*tmp(good,:)./h(good,ones(ncmp(i),1));
        
        % convert back
        if(amph(i))
            data(i).dep(:,1:2:end)=abs(tmp);
            data(i).dep(:,2:2:end)=angle(tmp);
        elseif(rlim(i))
            data(i).dep(:,1:2:end)=real(tmp);
            data(i).dep(:,2:2:end)=imag(tmp);
        else
            tmp=ifft(tmp,[],1,'symmetric');
            data(i).dep=tmp(1:npts(i),:);
        end
        
        % change class back
        data(i).dep=oclass(data(i).dep);
        
        % dep*
        depmen(i)=mean(data(i).dep(:)); 
        depmin(i)=min(data(i).dep(:)); 
        depmax(i)=max(data(i).dep(:));
    end
    
    % update header info
    data=changeheader(data,'scale',1,'idep',option.IDEP,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % since apply/remove sacpz bomb out so often...
    if(exist('i','var'))
        disp(sprintf('REMOVESACPZ bombed out on record: %d',i));
    end
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
