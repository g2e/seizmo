function [data,pz]=applysacpz(data,varargin)
%APPLYSACPZ    Applies SAC PoleZero responses to SEIZMO records
%
%    Usage:    [dataout,good]=applysacpz(datain)
%              [...]=applysacpz(...,'freqlimits',[f1 f2 f3 f4],...)
%              [...]=applysacpz(...,'units',UNITS,...)
%              [...]=applysacpz(...,'idep',IDEP,...)
%              [...]=applysacpz(...,'h2o',H2O,...)
%
%    Description: [DATAOUT,GOOD]=APPLYSACPZ(DATAIN) applies an instrument
%     response to records in DATAIN based on the associated SAC PoleZero
%     info.  If this info has not been added, then GETSACPZ is called.  If
%     any records in DATA do not have any associated SAC PoleZero info (as
%     is indicated by the .misc.has_sacpz struct field set by GETSACPZ),
%     the records are not returned in DATAOUT.  The secondary output, GOOD,
%     is a logical array indicating the records in DATAIN that had SAC
%     PoleZero info (.misc.has_sacpz set TRUE).  You may use a customized
%     PoleZero response on records by placing the info in the .misc.sacpz
%     struct field and making sure all records in DATA have the 
%     .misc.has_sacpz struct field set to TRUE or FALSE.  Otherwise
%     GETSACPZ is called on the entire dataset and the customized PoleZero
%     info is lost.  See GETSACPZ for info on the PoleZero layout.
%
%     [...]=APPLYSACPZ(...,'FREQLIMITS',[F1 F2 F3 F4],...) applies a
%     lowpass and a highpass taper that limits the spectrum of the
%     convolved records.  F1 and F2 give the highpass taper limits while
%     F3 and F4 specify the lowpass taper limits.  The highpass taper is
%     zero below F1 and unity above F2.  The lowpass taper is zero above F4
%     and unity below F3.  The tapers are cosine tapers applied in the
%     spectral domain.  This is an acausal filter and should not be used if
%     you want to preserve seismic phase onsets.
%
%     [...]=APPLYSACPZ(...,'UNITS',UNITS,...) sets the input units of
%     records.  By default, the SAC Polezero files will only give the
%     appropriate machine units when the response is applied if the input
%     records are in displacement in nanometers (a nanometers to meters
%     conversion is done internally before the response is applied).
%     Setting UNITS to 'vel' will properly convert any input records that
%     are in velocity to their instrument units based upon this assumption.
%     Valid strings are:
%      'none', 'dis', 'vel', and 'acc'.
%     The default units is 'none' (the same as 'dis').
%
%     [...]=APPLYSACPZ(...,'IDEP',IDEP,...) sets the output dependent
%     component label stored in the header field 'idep' to IDEP.  Typically
%     this is 'iunkn', 'ivolts', or 'icounts'.  IDEP may be a cell array of
%     strings to set each record separately.  The default is 'icounts'.
%
%     [...]=APPLYSACPZ(...,'H2O',H2O,...) applies a waterlevel factor to
%     the SAC PoleZero response.  H2O by default is 0, which has no effect.
%     This is provided for completeness with the option in REMOVESACPZ.
%
%    Notes:
%     - SAC PoleZero info should be set to convert machine units to
%       displacement in meters
%     - Output by default is machine units (IDEP is set to volts)
%     - the SCALE field is set to 1
%     - In order for GETSACPZ to identify the appropriate SAC PoleZero file
%       for each record, the following fields must be set correctly:
%        KNETWK, KSTNM, KHOLE, KCMPNM
%        NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC, B, E
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, IDEP, SCALE
%
%    Examples:
%     Apply instrument responses to velocity (nm/sec) records:
%      data=applysacpz(data,'units','vel');
%
%    See also: REMOVESACPZ, GETSACPZ, CONVOLVE

%     Version History:
%        Oct. 22, 2009 - initial version
%        Oct. 30, 2009 - added informative output on error
%        Feb.  3, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb. 16, 2010 - doc update, fixed related function list
%        May   5, 2010 - fixed upper frequency taper (thanks dsh)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   5, 2010 at 20:35 GMT

% todo:
% - update docs

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% import SEIZMO info
global SEIZMO

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
    % verbosity
    verbose=seizmoverbose;
    
    % attempt to access has_sacpz
    try
        pz=getsubfield([data.misc],'has_sacpz');
    catch
        % detail message
        if(verbose)
            disp('Not All Records Indicate PoleZero Status');
            disp('Attempting to Find SAC PoleZeros for All Records');
        end
        data=getsacpz(data);
        pz=getsubfield([data.misc],'has_sacpz');
    end
    
    % detail message
    if(verbose && any(~pz))
        disp(sprintf(['Record(s):\n' sprintf('%d ',find(~pz)) ...
            '\nDo Not Have SAC PoleZero Info.  Deleting!']));
    end
    
    % only use those with polezero info
    data=data(pz);
    
    % number of records
    nrecs=numel(data);
    
    % get header info
    leven=getlgc(data,'leven');
    [iftype,idep]=getenumid(data,'iftype','idep');
    [npts,ncmp,delta,e]=getheader(data,'npts','ncmp','delta','e');
    
    % get spectral
    rlim=strcmpi(iftype,'irlim');
    amph=strcmpi(iftype,'iamph');
    
    % get nyquist frequency
    nyq=1./(2*delta);
    if(any(rlim | amph)); nyq(rlim | amph)=e(rlim | amph); end
    
    % cannot do xyz records
    if(any(strcmpi(iftype,'ixyz')))
        error('seizmo:applysacpz:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(iftype,'ixyz')))
            '\nIllegal operation on XYZ record(s)!']);
    end
    
    % cannot do unevenly sampled records
    if(any(strcmpi(leven,'false')))
        error('seizmo:applysacpz:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false')))
            '\nInvalid operation on unevenly sampled records!']);
    end
    
    % valid values for strings
    % - this should be expanded to include all units
    % - need a matching wpow for each string
    % - and a matching idep unit
    valid.UNITS={...
        'n' 'none' ...
        'd' 'dis' 'disp' 'displacement' 'idisp' ...
        'v' 'vel' 'velo' 'velocity' 'ivel' ...
        'a' 'acc' 'accel' 'acceleration' 'iacc' ...
        'j' 'jerk' 'ijerk' ...
        's' 'snap' 'isnap' ...
        'c' 'crackle' 'icrackle' ...
        'p' 'pop' 'ipop' ...
        'absmnt' 'iabsmnt' 'absity' 'iabsity' 'abseler' 'iabseler' ...
        'abserk' 'iabserk' 'absnap' 'iabsnap' 'absackl' 'iabsackl' ...
        'abspop' 'iabspop' 'u' 'unkn' 'unknown' 'iunkn' 'volts' ...
        'ivolts' 'counts' 'icounts'};
    valid.WPOW=[0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 4 4 4 5 5 5 ...
        6 6 6 -1 -1 -2 -2 -3 -3 -4 -4 -5 -5 -6 -6 -7 -7 0 0 0 0 0 0 0 0];
    
    % get options from SEIZMO global
    ME=upper(mfilename);
    try
        varargin=[SEIZMO.(ME) varargin];
    catch
    end
    
    % default options
    flimbu=[-1*ones(nrecs,2) 2*nyq(:,[1 1])];
    varargin=[{'f' flimbu 'u' idep ...
        'id' 'icounts' 'h2o' zeros(nrecs,1)} varargin];
    
    % require all options to be strings
    if(~iscellstr(varargin(1:2:end)))
        error('seizmo:removesacpz:badInput',...
            'OPTIONS must be specified as strings!');
    end
    
    % check options
    for i=1:2:numel(varargin)
        value=varargin{i+1};
        
        % which option
        j=strmatch(lower(varargin{i}),{'freqlimits' 'units' 'idep' 'h2o'});
        switch j
            case 1 % freqlimits
                % assure real and correct size
                if(isreal(value) && any(size(value,1)==[1 nrecs]) ...
                        && size(value,2)<=4 && ndims(value)==2)
                    if(~isequal(value,sort(value,2)))
                        error('seizmo:applysacpz:badInput',...
                            ['FREQLIMITS must be [F1 F2 F3 F4]\n' ...
                            'where F1 <= F2 <= F3 <= F4!']);
                    end
                    if(size(value,1)==1)
                        flim=[value(ones(nrecs,1),:) ...
                            flimbu(:,(size(value,2)+1):4)];
                    else
                        flim=[value flimbu(:,(size(value,2)+1):4)];
                    end
                else
                    error('seizmo:applysacpz:badInput',...
                        'FREQLIMITS must be [F1 F2 F3 F4]!');
                end
            case 2 % units
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) ...
                        || ~any(numel(value)==[1 nrecs]) ...
                        || any(~ismember(value,valid.UNITS)))
                    error('seizmo:applysacpz:badInput',...
                        'UNITS must be ''DISP'' ''VEL'' or ''ACC''!');
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                units=value;
                
                % get associated WPOW
                [idx,idx]=ismember(units,valid.UNITS);
                wpow=valid.WPOW(idx);
            case 3 % idep
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:applysacpz:badInput',...
                        ['IDEP must be single string or have\n' ...
                        '1 string per record in DATA!']);
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                idep=value;
            case 4 % h2o
                if(~isreal(value) ...
                        || (~isscalar(value) && numel(value)~=nrecs) ...
                        || any(value<0))
                    error('seizmo:applysacpz:badInput',...
                        'H2O must be a real positive scalar or array!');
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                h2o=value;
            otherwise
                error('seizmo:applysacpz:badInput',...
                    'Unknown option: %s !',varargin{i});
        end
    end
    
    % detail message
    if(verbose)
        disp('Applying SAC PoleZero Response to Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    depmin=nan(nrecs,1); depmen=depmin; depmax=depmin;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % convert to complex spectra
        if(amph(i))
            nspts=npts(i);
            sdelta=delta(i);
            tmp=data(i).dep(:,1:2:end).*exp(1i*data(i).dep(:,2:2:end));
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
        taper2=taperfun('hann',nyq(i)-freq,nyq(i)-flim(i,[4 3])).';
        tmp=tmp.*taper1(:,ones(ncmp(1),1)).*taper2(:,ones(ncmp(1),1));
        
        % convert zpk to fap
        [a,p]=zpk2ap(freq,data(i).misc.sacpz.z,data(i).misc.sacpz.p,...
            data(i).misc.sacpz.k,wpow(i));
        h=((a+h2o(i)).*exp(1i*p)).';
        
        % apply response (over limited freqrange)
        % - divide by 1e9 to account for SAC PoleZero in meters
        tmp(good,:)=1e-9*tmp(good,:).*h(good,ones(ncmp(i),1));
        
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
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update header info
    data=changeheader(data,'scale',1,'idep',idep,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % since apply/remove sacpz bomb out so often...
    if(exist('i','var'))
        disp(sprintf('APPLYSACPZ bombed out on record: %d',i));
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
