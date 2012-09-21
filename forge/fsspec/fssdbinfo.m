function [maxdb,meddb,mindb]=fssdbinfo(s,fflag,varargin)
%FSSDBINFO    Returns the min/median/max dB for a fss struct
%
%    Usage:    [maxdb,meddb,mindb]=fssdbinfo(s)
%              [...]=fssdbinfo(s,fflag)
%              [...]=fssdbinfo(s,fflag,frng,bazrng,srng,esrng,nsrng)
%
%    Description:
%     [MAXDB,MEDDB,MINDB]=FSSDBINFO(S) returns the decibel limits and
%     median of the frequency-slowness spectra(s) in the fss struct S.
%     This is useful for quick identification for strong plane wave
%     coherency.  The outputs are structs with the following format:
%       .db    == decibel value
%       .slow  == magnitude of the horizontal slowness (in sec/deg)
%       .baz   == backazimuth (in degrees)
%       .freq  == frequency (in Hz)
%     Each field is equal in size to the input struct S.  So MAXDB.db(3),
%     MAXDB.slow(3), MAXDB.baz(3), & MAXDB.freq(3) give info about the max
%     peak in S(3).  Please note that MEDDB does not return any location
%     info as a median's location is not useful/straightfoward.
%
%     [...]=FSSDBINFO(S,FFLAG) sets if each frequency is processed
%     individually.  If FFLAG is TRUE then each frequency is processed
%     separately and the output will have NFREQ entries.  The default is
%     FALSE.
%
%     [...]=FSSDBINFO(S,FFLAG,FRNG,BAZRNG,SRNG,ESRNG,NSRNG) returns the dB
%     info about a particular subsection of the frequency-slowness
%     spectra(s).  FRNG is the frequency range in Hz.  BAZRNG is the back-
%     azimuth range in degrees.  SRNG is the horizontal slowness magnitude
%     range in sec/deg.  ESRNG & NSRNG are the East & North slowness ranges
%     in sec/deg.  The defaults for all ranges do not exclude any portion
%     of the spectra.
%
%    Notes:
%
%    Examples:
%     % Get a dispersion curve for Rayleigh wave noise:
%     s=fss(data,50,301,[1/100 1/10],true);
%     db=fssdbinfo(s,true,[],[25 35]);
%     figure; plot(db.freq,db.slow);
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, FSSSUB, FSSAVG

%     Version History:
%        May  26, 2010 - initial version
%        June 30, 2010 - added range arguments
%        July  1, 2010 - range bugfix
%        July  6, 2010 - update for new struct
%        July 12, 2010 - baz range issue bugfix (for sane ranges)
%        July 16, 2010 - output structs with db point info
%        Mar. 29, 2012 - minor doc update
%        Sep. 13, 2012 - adapted from fkdbinfo & geofssdbinfo
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2012 at 14:25 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check fss struct
error(chkfss(s));

% default flag
if(nargin<2 || isempty(fflag)); fflag=false; end

% subspectra if frng passed
if(nargin>2 && ~isempty(varargin{1}))
    s=fsssub(s,varargin{1});
end

% loop over every element
sz=size(s);
mindb=struct('db',[],'slow',[],'baz',[],'freq',[]);
mindb=repmat(mindb,sz); meddb=mindb; maxdb=mindb;
for i=1:numel(s)
    % get grids
    [baz,smag]=fssgrids(s(i));
    sz=size(s(i).spectra);
    if(numel(sz)==2); sz(3)=1; end
    
    % scale factor if unwhitened
    if(~s(i).whiten); s(i).spectra=s(i).spectra/(2*pi); end
    
    % convert to dB
    switch s(i).method
        case {'coarray' 'full'}
            s(i).spectra=10*log10(abs(real(s(i).spectra)));
        otherwise
            s(i).spectra=10*log10(s(i).spectra);
    end
    
    % force freq field to column vector
    s(i).freq=s(i).freq(:);
    
    % stats
    if(fflag)
        % reduce spectra
        if(sz(3)~=1)
            % need map indexing so squash frequency dimension
            idx=fssidx(fssavg(s),[],varargin{2:end});
        else
            idx=fssidx(s,[],varargin{2:end});
        end
        s(i).spectra=permute(s(i).spectra,[3 1 2]);
        s(i).spectra=s(i).spectra(:,idx(:));
        idx=find(idx);
        
        % get values for each frequency in that spectra
        meddb(i).db=nanmedian(s(i).spectra,2);
        meddb(i).freq=s(i).freq;
        [mindb(i).db,midx]=min(s(i).spectra,[],2);
        mindb(i).slow=smag(idx(midx));
        mindb(i).baz=baz(idx(midx));
        mindb(i).freq=s(i).freq;
        [maxdb(i).db,midx]=max(s(i).spectra,[],2);
        maxdb(i).slow=smag(idx(midx));
        maxdb(i).baz=baz(idx(midx));
        maxdb(i).freq=s(i).freq;
    else
        % reduce spectra
        idx=fssidx(s,[],varargin{2:end});
        [r,c,p]=ind2sub(sz,find(idx));
        s(i).spectra=s(i).spectra(idx);
        
        % get values for that spectra
        meddb(i).db=nanmedian(s(i).spectra);
        [mindb(i).db,idx]=min(s(i).spectra);
        mindb(i).slow=smag(r(idx),c(idx));
        mindb(i).baz=baz(r(idx),c(idx));
        mindb(i).freq=s(i).freq(p(idx));
        [maxdb(i).db,idx]=max(s(i).spectra);
        maxdb(i).slow=smag(r(idx),c(idx));
        maxdb(i).baz=baz(r(idx),c(idx));
        maxdb(i).freq=s(i).freq(p(idx));
    end
    
    % deal with freq averaging
    if(sz(3)==1)
        meddb(i).freq=s(i).freq.';
        mindb(i).freq=s(i).freq.';
        maxdb(i).freq=s(i).freq.';
    end
end

end
