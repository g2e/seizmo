function [data,pz]=applysacpz(data,varargin)
%APPLYSACPZ    Applies SAC PoleZero responses to SEIZMO records
%
%    Usage:    [dataout,good]=applysacpz(datain)
%              [...]=applysacpz(...,'taperlimits',[f1 f2 f3 f4],...)
%              [...]=applysacpz(...,'tapertype',TYPE,...)
%              [...]=applysacpz(...,'taperopt',VALUE,...)
%              [...]=applysacpz(...,'units',UNITS,...)
%              [...]=applysacpz(...,'idep',IDEP,...)
%              [...]=applysacpz(...,'h2o',H2O,...)
%
%    Description:
%     [DATAOUT,GOOD]=APPLYSACPZ(DATAIN) applies an instrument response to
%     records in DATAIN based on the associated SAC PoleZero info.  If this
%     info has not been added, then GETSACPZ is called.  If any records in
%     DATA do not have any associated SAC PoleZero info (as is indicated by
%     the .misc.has_sacpz struct field set by GETSACPZ), the records are
%     not returned in DATAOUT.  The secondary output, GOOD, is a logical
%     array indicating the records in DATAIN that had SAC PoleZero info
%     (.misc.has_sacpz set TRUE).  You may use a customized PoleZero
%     response on records by placing the info in the .misc.sacpz struct
%     field and making sure all records in DATA have the .misc.has_sacpz
%     struct field set to TRUE or FALSE.  Otherwise GETSACPZ is called on
%     the entire dataset and the customized PoleZero info is lost.  See
%     GETSACPZ for info on the PoleZero layout.  Please note that the
%     polezero response is expected to convert records in displacement
%     (meters) to machine units (counts).  The response is adjusted for
%     records that are in ground units other than displacement (as
%     indicated by the IDEP header field) so they are in the same machine
%     units (counts) as a displacement would be.  Also note that records
%     are assumed to be nanometers-based (which does not match the polezero
%     assumption) and so they are scaled internally to correct for this.
%
%     [...]=APPLYSACPZ(...,'TAPERLIMITS',[F1 F2 F3 F4],...) applies a
%     lowpass and a highpass taper that limits the spectrum of the
%     convolved records.  F1 and F2 give the highpass taper limits while
%     F3 and F4 specify the lowpass taper limits.  The highpass taper is
%     zero below F1 and unity above F2.  The lowpass taper is zero above F4
%     and unity below F3.  The tapers are cosine tapers applied in the
%     spectral domain.  This is an acausal filter and should not be used if
%     you want to preserve seismic phase onsets.
%
%     [...]=APPLYSACPZ(...,'TAPERTYPE',TYPE,...) alters the taper type of
%     the frequency limiting.  The default type is 'hann' and is a cosine
%     taper.  TYPE may be any offered by the TAPERFUN function.
%
%     [...]=APPLYSACPZ(...,'TAPEROPT',VALUE,...) adjusts the taper option
%     if there is one.  See TAPERFUN for details.
%
%     [...]=APPLYSACPZ(...,'UNITS',UNITS,...) overrides the units found in
%     the IDEP header field with UNITS.  This could also be done by setting
%     the IDEP field for the records prior to calling APPLYSACPZ.  This is
%     mainly useful for handling special cases like polezero info that
%     is not set to convert between displacement & counts (see the Examples
%     section below).  The default is found from the records' IDEP field.
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
%     - SAC PoleZero info should be set to convert between machine units &
%       displacement in meters
%     - Records should be in nanometers-based ground units (records in
%       non-ground units are treated as displacement records in nanometers) 
%     - Output by default is machine units (IDEP is set to counts)
%     - the SCALE field is set to 1
%     - In order for GETSACPZ to identify the appropriate SAC PoleZero file
%       for each record, the following fields must be set correctly:
%        KNETWK, KSTNM, KHOLE, KCMPNM
%        NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC, B, E
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, IDEP, SCALE
%
%    Examples:
%     % Apply instrument responses to records in any known ground unit:
%     data=applysacpz(data);
%
%     % Say your response info converts between velocity & counts.  Assume
%     % your records are in velocity and you want them in counts.  Using
%     % the default will improperly convert your records to displacement as
%     % part of the response removal (because we assume the response
%     % converts between displacement & counts).  To get the proper action
%     % set UNITS to 'dis' so no adjustment is made:
%     data=applysacpz(data,'units','dis');
%
%     % These should all give nearly exactly the same result (error is due
%     % to the discrete integration/differentiation):
%     data=applysacpz(data);
%     data=applysacpz(differentiate(data));
%     data=applysacpz(integrate(data));
%
%     % My records are in meters!  What do I do?  Convert to nanometers:
%     data=applysacpz(divide(data,1e9));
%
%    See also: REMOVESACPZ, GETSACPZ, CONVOLVE

%     Version History:
%        Oct. 22, 2009 - initial version
%        Oct. 30, 2009 - added informative output on error
%        Feb.  3, 2010 - seizmoverbose support, proper SEIZMO handling
%        Feb. 16, 2010 - doc update, fixed related function list
%        May   5, 2010 - fixed upper frequency taper (thanks dsh)
%        May   7, 2010 - doc update, changed global option passing,
%                        can now pass partial option strings, fix bug in
%                        idep/units, allow many more ground units
%        Aug. 19, 2010 - removed ifft symmetric flag and do real conversion
%                        afterwards, no longer use strmatch for options
%        Aug. 20, 2010 - taperopt/tapertype options added, better
%                        checkheader usage
%        Aug. 25, 2010 - drop SEIZMO global, fix taper option bug
%        June  9, 2011 - changing freqlimits to taperlimits, output taper
%                        limits to terminal, fixed bomb-out bug
%        Feb.  3, 2012 - doc update
%        Mar. 13, 2012 - use getheader improvements
%        May  30, 2012 - pow2pad=0 by default
%        Mar. 10, 2014 - works with new sacpz format, error if no sacpz,
%                        reduced computations (skip neg freq), skip/delete
%                        records with bad responses, fixed upper taper eps
%                        bug that could cause flatlining
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 10, 2014 at 20:30 GMT

% todo:
% - standard responses
% - maybe we should just have a wpow option rather than units
% - meters/nanometers flag
% - scale field should be ignored rather than changed

% check nargin
error(nargchk(1,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data,...
        'FALSE_LEVEN','ERROR',...
        'XYZ_IFTYPE','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
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
        fprintf(['Record(s):\n' sprintf('%d ',find(~pz)) ...
            '\nDo Not Have SAC PoleZero Info.  Deleting!\n']);
    end
    
    % only use those with polezero info
    data=data(pz);
    
    % number of records
    nrecs=numel(data);
    
    % handle no records with pz
    if(nrecs==0)
        error('seizmo:applysacpz:noPoleZeros',...
            'No SAC PoleZeros are available for records in DATA!');
    end
    
    % now remove those with a bad sacpz
    pz=find(pz);
    keep=true(nrecs,1);
    for i=1:nrecs
        if(isfield(data(i).misc.sacpz,'bad') && data(i).misc.sacpz.bad)
            warning('seizmo:removesacpz:badSACPZ',...
                'Record %d has a bad SAC PoleZero response! Deleting!',...
                pz(i));
            keep(i)=false;
        end
    end
    data=data(keep);
    
    % number of records
    nrecs=numel(data);
    
    % handle no records with pz
    if(nrecs==0)
        error('seizmo:removesacpz:noPoleZeros',...
            'No Good SAC PoleZeros are available for records in DATA!');
    end
    
    % get header info
    [npts,ncmp,delta,e,iftype,idep]=getheader(data,...
        'npts','ncmp','delta','e','iftype id','idep id');
    
    % get spectral
    rlim=strcmpi(iftype,'irlim');
    amph=strcmpi(iftype,'iamph');
    
    % get nyquist frequency
    nyq=1./(2*delta);
    if(any(rlim | amph)); nyq(rlim | amph)=e(rlim | amph); end
    
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
    
    % default options
    tlimbu=[-1*ones(nrecs,2) 2*nyq 2*nyq+10*eps(2*nyq)]; tlim=tlimbu;
    tlimstr={'-1' '-1' '2*NYQUIST' '2*NYQUIST'};
    varargin=[{'t' 'hann' 'o' [] 'u' idep ...
        'id' 'icounts' 'h2o' zeros(nrecs,1)} varargin];
    
    % require all options to be strings
    if(~iscellstr(varargin(1:2:end)))
        error('seizmo:removesacpz:badInput',...
            'OPTIONS must be specified as strings!');
    end
    
    % check options
    for z=1:2:numel(varargin)
        value=varargin{z+1};
        
        % which option
        switch lower(varargin{z})
            case {'fl' 'freq' 'fr' 'freqlim' 'freql' 'freqlimits' 'frq' ...
                    'freqlimit' 'frqlim' 'f' 'tl' 'ta' 'tap' 'tpr' ...
                    'taperlim' 'taperl' 'taperlimits' 'taperlimit' ...
                    'tprlim' 'l'}
                % assure real and correct size
                if(isreal(value) && any(size(value,1)==[1 nrecs]) ...
                        && size(value,2)<=4 && ndims(value)==2)
                    % make sure it is sorted
                    if(~isequal(value,sort(value,2)))
                        error('seizmo:applysacpz:badInput',...
                            ['TAPERLIMITS must be [F1 F2 F3 F4]\n' ...
                            'where F1 <= F2 <= F3 <= F4!']);
                    end
                    
                    % expand to Nx4 as needed
                    if(size(value,1)==1)
                        tlim=[value(ones(nrecs,1),:) ...
                            tlimbu(:,(size(value,2)+1):4)];
                    else
                        tlim=[value tlimbu(:,(size(value,2)+1):4)];
                    end
                    
                    % update taperlimits strings
                    if(size(value,1)==1)
                        for a=1:numel(value)
                            tlimstr{a}=num2str(value(a),'%g');
                        end
                    else % check for variable limits (use 'VARIABLE')
                        for a=1:size(value,2)
                            if(isscalar(unique(value(:,a))))
                                tlimstr{a}=num2str(value(1,a),'%g');
                            else
                                tlimstr{a}='VARIABLE';
                            end
                        end
                    end
                else
                    error('seizmo:applysacpz:badInput',...
                        'TAPERLIMITS must be [F1 F2 F3 F4]!');
                end
            case {'u' 'un' 'unit' 'units'}
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
            case {'i' 'id' 'dep' 'idep'}
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:applysacpz:badInput',...
                        ['IDEP must be single string or have\n' ...
                        '1 string per record in DATA!']);
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                idep=value;
            case {'h' 'h2' 'h2o'}
                if(~isreal(value) ...
                        || (~isscalar(value) && numel(value)~=nrecs) ...
                        || any(value<0))
                    error('seizmo:applysacpz:badInput',...
                        'H2O must be a real positive scalar or array!');
                end
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                h2o=value;
            case {'tapertype' 't' 'tt' 'ttype'}
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:applysacpz:badInput',...
                        ['TAPERTYPE must be a single string or have\n' ...
                        '1 string per record in DATA!']);
                end
                ttype=value;
                if(isscalar(ttype)); ttype(1:nrecs,1)=ttype; end
            case {'taperoption' 'o' 'to' 'topt'}
                if(isempty(value)); topt=cell(nrecs,1); continue; end
                if(~isreal(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:applysacpz:badInput',...
                        'TAPEROPT must be a real-valued scalar or array!');
                end
                topt=value;
                if(isscalar(topt)); topt(1:nrecs,1)=topt; end
            otherwise
                error('seizmo:applysacpz:badInput',...
                    'Unknown option: %s !',varargin{z});
        end
    end
    
    % detail message
    if(verbose)
        disp('Applying SAC PoleZero Response to Record(s)');
        fprintf('Taper FullStop Limits: %10s Hz, %10s Hz\n',...
            tlimstr{1},tlimstr{4})
        fprintf('Taper FullPass Limits: %10s Hz, %10s Hz\n',...
            tlimstr{2},tlimstr{3})
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
            tmp=data(i).dep(:,1:2:end).*exp(1i*data(i).dep(:,2:2:end));
        elseif(rlim(i))
            nspts=npts(i);
            tmp=complex(data(i).dep(:,1:2:end),data(i).dep(:,2:2:end));
        else
            nspts=2^nextpow2(npts(i));
            tmp=fft(data(i).dep,nspts,1); % no need to scale by delta
        end
        
        % just the positive frequencies
        freq=linspace(0,nyq(i),nspts/2+1);
        
        % convert zpk to complex (this is several orders better than SAC)
        h=zpk2cmplx(freq,data(i).misc.sacpz.z{1},...
            data(i).misc.sacpz.p{1},data(i).misc.sacpz.k,wpow(i));
        h=((abs(h)+h2o(i)).*exp(1i*angle(h))).'; % add waterlevel...
        
        % now apply taper
        h=h.*taperfun(ttype{i},freq,tlim(i,1:2),topt(i)).'...
            .*taperfun(ttype{i},freq,tlim(i,[4 3]),topt(i)).';
        h(isinf(h) | isnan(h))=0; % care for trouble spots in response
        
        % indexing to only work within the taper limits
        good=freq>=tlim(i,1) & freq<=tlim(i,4);
        tmpgood=[good false(1,nspts/2-1)];
        
        % apply transfer function
        % - divide by 1e9 to account for SAC PoleZero in meters
        tmp(tmpgood,:)=1e-9*tmp(tmpgood,:).*h(good,ones(ncmp(i),1));
        tmp(~tmpgood,:)=0; % no response outside taper limits
        tmp(nspts/2+2:end,:)=conj(tmp(nspts/2:-1:2,:)); % neg frequencies
        tmp(1)=0; % force 0Hz to 0
        tmp(nspts/2+1)=abs(tmp(nspts/2+1)); % no constraint on NyqHz phase
        %tmp(nspts)=abs(tmp(nspts)); % SAC transfer bug
        
        % convert back
        if(amph(i))
            data(i).dep(:,1:2:end)=abs(tmp);
            data(i).dep(:,2:2:end)=angle(tmp);
        elseif(rlim(i))
            data(i).dep(:,1:2:end)=real(tmp);
            data(i).dep(:,2:2:end)=imag(tmp);
        else
            tmp=real(ifft(tmp,[],1));
            data(i).dep=tmp(1:npts(i),:);
        end
        
        % change class back
        data(i).dep=oclass(data(i).dep);
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:)); 
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
        fprintf('APPLYSACPZ bombed out on record: %d\n',i);
        data(i)
        data(i).misc.sacpz
        data(i).misc.sacpz.z{1}
        data(i).misc.sacpz.p{1}
        data(i).misc.sacpz.k
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
