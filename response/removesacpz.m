function [data,pz]=removesacpz(data,varargin)
%REMOVESACPZ    Removes SAC PoleZero responses from SEIZMO records
%
%    Usage:    [dataout,good]=removesacpz(datain)
%              [...]=removesacpz(...,'taperlimits',[f1 f2 f3 f4],...)
%              [...]=removesacpz(...,'tapertype',TYPE,...)
%              [...]=removesacpz(...,'taperopt',VALUE,...)
%              [...]=removesacpz(...,'units',UNITS,...)
%              [...]=removesacpz(...,'idep',IDEP,...)
%              [...]=removesacpz(...,'h2o',H2O,...)
%
%    Description:
%     [DATAOUT,GOOD]=REMOVESACPZ(DATAIN) removes the instrument response
%     from records in DATAIN based on the associated SAC PoleZero info.  If
%     this info has not been added, then GETSACPZ is called.  If any
%     records in DATA do not have any associated SAC PoleZero info (as is
%     set in the .misc.has_sacpz struct field by GETSACPZ), the records are
%     not returned in DATAOUT.  The secondary output, GOOD, is a logical
%     array indicating the records in DATAIN that had SAC PoleZero info
%     (.misc.has_sacpz set TRUE).  You may use a customized PoleZero
%     response on records by placing the info in the .misc.sacpz struct
%     field and making sure _ALL_ records in DATA have the .misc.has_sacpz
%     struct field set to TRUE or FALSE.  See GETSACPZ for info on the
%     PoleZero layout.
%
%     [...]=REMOVESACPZ(...,'TAPERLIMITS',[F1 F2 F3 F4],...) applies a
%     lowpass and a highpass taper that limits the spectrum of the
%     deconvolved records.  F1 and F2 give the highpass taper limits while
%     F3 and F4 specify the lowpass taper limits.  The highpass taper is
%     zero below F1 and unity above F2.  The lowpass taper is zero above F4
%     and unity below F3.  The tapers are cosine tapers applied in the
%     spectral domain.  This is an acausal filter and should not be used if
%     you want to preserve seismic phase onsets.  Use option H2O to
%     stabilize the deconvolution without frequency domain tapering.
%     Default TAPERLIMITS are [-1 -1 2*NYQ 2*NYQ] where NYQ is the nyquist
%     frequency for each particular record.  The defaults will not apply
%     any tapering to the records.  Note that TAPERLIMITS may be specified
%     also as [F1 F2] or [F1 F2 F3].  In these cases the remaining unset
%     values are kept at their defaults.
%
%     [...]=REMOVESACPZ(...,'TAPERTYPE',TYPE,...) alters the taper type of
%     the frequency limiting.  The default type is 'hann' and is a cosine
%     taper.  TYPE may be any offered by the TAPERFUN function.
%
%     [...]=REMOVESACPZ(...,'TAPEROPT',VALUE,...) adjusts the taper option
%     if there is one.  See TAPERFUN for details.
%
%     [...]=REMOVESACPZ(...,'UNITS',UNITS,...) allows changing the specific
%     ground units that the records are converted to (displacement in nm,
%     velocity in nm/sec, or acceleration in nm/sec^2).  UNITS must be one
%     of the following strings: 'DISP', 'VEL', or 'ACC'.  Note that this
%     will automatically update the IDEP values to match (next option).
%     Default UNITS is 'DISP'.
%
%     [...]=REMOVESACPZ(...,'IDEP',IDEP,...) sets the output dependent
%     component label stored in the header field 'idep' to IDEP.  Typically
%     this is 'idisp', 'ivel', 'iacc', or 'iunkn'.  IDEP may be a cell
%     array of strings to set each record separately.  The default IDEP is
%     'idisp'.  Note that using the 'UNITS' option after the 'IDEP' option
%     will replace the values set by IDEP with those corresponding to
%     UNITS.
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
%     % Remove the instrument response for a dataset, converting to
%     % velocity and ignoring all periods at >200s while tapering periods
%     % from 200 to 150 seconds:
%     data=removesacpz(data,'units','vel','taperlimits',[1/200 1/150]);
%
%    See also: APPLYSACPZ, GETSACPZ, DECONVOLVE, ZPK2CMPLX, ZPK2AP

%     Version History:
%        Oct. 22, 2009 - initial version
%        Oct. 30, 2009 - added informative output on error
%        Feb.  3, 2010 - proper SEIZMO handling
%        Feb. 16, 2010 - added info about custom PoleZero usage
%        May   5, 2010 - cleaned up documentation, fixed upper frequency
%                        taper (thanks dsh), changed global option passing,
%                        can now pass partial option strings, fix bug in
%                        idep/units, allow many more ground units, handle
%                        0 response returning nans (set nans to 0)
%        Aug. 19, 2010 - removed ifft symmetric flag and do real conversion
%                        afterwards, no longer use strmatch for options
%        Aug. 20, 2010 - taperopt/tapertype options added, better
%                        checkheader usage
%        Aug. 25, 2010 - drop SEIZMO global, fix taper option bug
%        Sep. 28, 2010 - more info when response bombs out
%        June  9, 2011 - changing freqlimits to taperlimits, output taper
%                        limits to terminal, fixed bomb-out bug
%        Feb.  3, 2012 - doc update
%        Mar. 13, 2012 - use getheader improvements
%        May  30, 2012 - pow2pad=0 by default
%        June  1, 2012 - reduced computions (skip neg freq), found sac bug
%        Mar. 10, 2014 - works with new sacpz format, error if no sacpz,
%                        skip/delete records with bad responses, doc
%                        update, fixed upper taper eps bug that could cause
%                        flatlining
%        Apr. 28, 2014 - renamed loop increment variable for debugging
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 28, 2014 at 20:30 GMT

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
            disp('Not All Records Indicate PoleZero Status.');
            disp('Attempting to Find SAC PoleZeros for All Records.');
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
        error('seizmo:removesacpz:noPoleZeros',...
            'No SAC PoleZeros are available for records in DATA!');
    end
    
    % now remove those with a bad sacpz
    pz=find(pz);
    keep=true(nrecs,1);
    for j=1:nrecs
        if(isfield(data(j).misc.sacpz,'bad') && data(j).misc.sacpz.bad)
            warning('seizmo:removesacpz:badSACPZ',...
                'Record %d has a bad SAC PoleZero response! Deleting!',...
                pz(j));
            keep(j)=false;
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
    [npts,ncmp,delta,e,iftype]=getheader(data,...
        'npts','ncmp','delta','e','iftype id');
    
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
    valid.IDEP={'idisp' 'idisp' 'idisp' 'idisp' 'idisp' 'idisp' 'idisp' ...
        'ivel' 'ivel' 'ivel' 'ivel' 'ivel' 'iacc' 'iacc' 'iacc' 'iacc' ...
        'iacc' 'ijerk' 'ijerk' 'ijerk' 'isnap' 'isnap' 'isnap' ...
        'icrackle' 'icrackle' 'icrackle' 'ipop' 'ipop' 'ipop' 'iabsmnt' ...
        'iabsmnt' 'iabsity' 'iabsity' 'iabseler' 'iabseler' 'iabserk' ...
        'iabserk' 'iabsnap' 'iabsnap' 'iabsackl' 'iabsackl' 'iabspop' ...
        'iabspop' 'iunkn' 'iunkn' 'iunkn' 'iunkn' 'ivolts' 'ivolts' ...
        'icounts' 'icounts'};
    
    % default options
    tlimbu=[-1*ones(nrecs,2) 2*nyq 2*nyq+10*eps(2*nyq)]; tlim=tlimbu;
    tlimstr={'-1' '-1' '2*NYQUIST' '2*NYQUIST'};
    varargin=[{'t' 'hann' 'o' [] 'u' 'dis' ...
        'id' 'idisp' 'h2o' zeros(nrecs,1)} varargin];
    
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
                        error('seizmo:removesacpz:badInput',...
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
                    error('seizmo:removesacpz:badInput',...
                        'TAPERLIMITS must be [F1 F2 F3 F4]!');
                end
            case {'u' 'un' 'unit' 'units'}
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) ...
                        || ~any(numel(value)==[1 nrecs]) ...
                        || any(~ismember(value,valid.UNITS)))
                    error('seizmo:removesacpz:badInput',...
                        'UNITS must be ''DISP'' ''VEL'' or ''ACC''!');
                end
                
                % expand scalars
                if(isscalar(value)); value=value(ones(nrecs,1),1); end
                units=value;
                
                % get associated WPOW/IDEP
                [idx,idx]=ismember(units,valid.UNITS);
                wpow=valid.WPOW(idx);
                idep=valid.IDEP(idx);
            case {'i' 'id' 'dep' 'idep'}
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:removesacpz:badInput',...
                        ['IDEP must be a single string or have\n' ...
                        '1 string per record in DATA!']);
                end
                idep=value;
                if(isscalar(idep)); idep(1:nrecs,1)=idep; end
            case {'h' 'h2' 'h2o'}
                if(~isreal(value) || ~any(numel(value)==[1 nrecs]) ...
                        || any(value<0))
                    error('seizmo:removesacpz:badInput',...
                        'H2O must be a real positive scalar or array!');
                end
                h2o=value;
                if(isscalar(h2o)); h2o(1:nrecs,1)=h2o; end
            case {'tapertype' 't' 'tt' 'ttype'}
                if(ischar(value)); value=cellstr(value); end
                if(~iscellstr(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:removesacpz:badInput',...
                        ['TAPERTYPE must be a single string or have\n' ...
                        '1 string per record in DATA!']);
                end
                ttype=value;
                if(isscalar(ttype)); ttype(1:nrecs,1)=ttype; end
            case {'taperoption' 'o' 'to' 'topt'}
                if(isempty(value)); topt=cell(nrecs,1); continue; end
                if(~isreal(value) || ~any(numel(value)==[1 nrecs]))
                    error('seizmo:removesacpz:badInput',...
                        'TAPEROPT must be a real-valued scalar or array!');
                end
                topt=value;
                if(isscalar(topt)); topt(1:nrecs,1)=topt; end
            otherwise
                error('seizmo:removesacpz:badInput',...
                    'Unknown option: %s !',varargin{z});
        end
    end
    
    % detail message
    if(verbose)
        disp('Removing SAC PoleZero Response from Record(s)');
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
        
        % convert data to complex spectra
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
        h=1./h.*taperfun(ttype{i},freq,tlim(i,1:2),topt(i)).'...
            .*taperfun(ttype{i},freq,tlim(i,[4 3]),topt(i)).';
        h(isinf(h) | isnan(h))=0; % care for trouble spots in response
        
        % indexing to only work within the taper limits
        good=freq>=tlim(i,1) & freq<=tlim(i,4);
        tmpgood=[good false(1,nspts/2-1)];
        
        % apply transfer function
        % - multiply by 1e9 to account for SAC PoleZero in meters
        tmp(tmpgood,:)=1e9*tmp(tmpgood,:).*h(good,ones(ncmp(i),1));
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
        fprintf('REMOVESACPZ bombed out on record: %d\n',i);
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
