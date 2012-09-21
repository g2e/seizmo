function [maxdb,meddb,mindb]=geofssdbinfo(s,fflag,sflag)
%GEOFSSDBINFO    Returns the min/median/max dB for a freq-slow spectra
%
%    Usage:    [maxdb,meddb,mindb]=geofssdbinfo(s)
%              [maxdb,meddb,mindb]=geofssdbinfo(s,fflag,sflag)
%
%    Description:
%     [MAXDB,MEDDB,MINDB]=GEOFSSDBINFO(S) returns the decibel limits and
%     median of the frequency-slowness-position spectra(s) in the GEOFSS
%     struct S.  This is useful for identifying waves with strong spherical
%     wave coherency.  The outputs are structs with the following format:
%       .db        == decibel value
%       .slow      == magnitude of the horizontal slowness (in sec/deg)
%       .latlon    == [latitude longitude] (in degrees)
%       .freq      == frequency (in Hz)
%     The output is equal in size to the input struct S.  So MAXDB(3).db,
%     MAXDB(3).slow, MAXDB(3).latlon, & MAXDB(3).freq give info about the
%     max peak in S(3).  Please note that MEDDB does not return any
%     location info as a median's location is not useful/straightfoward.
%
%     [MAXDB,MEDDB,MINDB]=GEOFSSDBINFO(S,FFLAG,SFLAG) sets if each
%     frequency and each slowness is processed individually.  If FFLAG is
%     TRUE then each frequency is processed separately and the output will
%     have NFREQ entries.  If SFLAG is TRUE then each slowness is processed
%     separately and the output will have NSLOW entries.  If both are TRUE
%     then the outputs will have NFREQxNSLOW entries.
%
%    Notes:
%
%    Examples:
%     % A quick way to look at how slowness varies with frequency:
%     db=geofssdbinfo(s,true);
%     plot(db.freq,db.slow);
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSSUB, GEOFSSAVG

%     Version History:
%        July 12, 2010 - initial version
%        July 18, 2010 - output structs with db point info
%        Mar. 24, 2012 - minor doc update
%        June 11, 2012 - adapted from geofkdbinfo, drop sub inputs, add
%                        flag inputs, output changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2012 at 22:25 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check struct
error(chkgeofss(s));

% default flags
if(nargin<2 || isempty(fflag)); fflag=false; end
if(nargin<3 || isempty(sflag)); sflag=false; end

% check flags
if(~isscalar(fflag) || ~islogical(fflag))
    error('seizmo:geofssdbinfo:badInput',...
        'FFLAG must be TRUE or FALSE!');
elseif(~isscalar(sflag) || ~islogical(sflag))
    error('seizmo:geofssdbinfo:badInput',...
        'SFLAG must be TRUE or FALSE!');
end

% flags vs vector
v=reshape([s.vector],[],2);
if(fflag && any(~v(:,1)))
    error('seizmo:geofssdbinfo:badInput',...
        'FFLAG must be FALSE for frequency averaged spectra!');
elseif(sflag && any(~v(:,2)))
    error('seizmo:geofssdbinfo:badInput',...
        'SFLAG must be FALSE for slowness averaged spectra!');
end

% loop over every element
sz=size(s);
mindb=struct('db',[],'slow',[],'latlon',[],'freq',[]);
mindb=repmat(mindb,sz); meddb=mindb; maxdb=mindb;
for i=1:numel(s)
    % convert to dB
    switch s(i).method
        case {'coarray' 'full'}
            s(i).spectra=10*log10(abs(real(s(i).spectra)));
        otherwise
            s(i).spectra=10*log10(s(i).spectra);
    end
    
    % force freq field to column vector
    s(i).freq=s(i).freq(:);
    
    % handle based on flags
    if(fflag && sflag)
        % get values for each freq & slow
        sz=size(s(i).spectra);
        meddb(i).db=median(s(i).spectra(:,:)).';
        [mindb(i).db,idx]=min(s(i).spectra(:,:));
        mindb(i).db=mindb(i).db.';
        mindb(i).latlon=s(i).latlon(idx,:);
        mindb(i).slow=repmat(s(i).slow,sz(3),1);
        mindb(i).freq=repmat(s(i).freq,1,sz(2)).';
        mindb(i).freq=mindb(i).freq(:);
        [maxdb(i).db,idx]=max(s(i).spectra(:,:));
        maxdb(i).db=maxdb(i).db.';
        maxdb(i).latlon=s(i).latlon(idx,:);
        maxdb(i).slow=repmat(s(i).slow,sz(3),1);
        maxdb(i).freq=repmat(s(i).freq,1,sz(2)).';
        maxdb(i).freq=mindb(i).freq(:);
    elseif(fflag)
        % get values for each frequency
        sz=size(s(i).spectra);
        s(i).spectra=permute(s(i).spectra,[3 1 2]);
        meddb(i).db=median(s(i).spectra(:,:),2);
        [mindb(i).db,idx]=min(s(i).spectra(:,:),[],2);
        [r,c]=ind2sub(sz(1:2),idx);
        mindb(i).latlon=s(i).latlon(r,:);
        mindb(i).slow=s(i).slow(c);
        mindb(i).freq=s(i).freq;
        [maxdb(i).db,idx]=max(s(i).spectra(:,:),[],2);
        [r,c]=ind2sub(sz(1:2),idx);
        maxdb(i).latlon=s(i).latlon(r,:);
        maxdb(i).slow=s(i).slow(c);
        maxdb(i).freq=s(i).freq;
    elseif(sflag)
        % get values for each slowness
        sz=size(s(i).spectra);
        s(i).spectra=permute(s(i).spectra,[2 1 3]);
        meddb(i).db=median(s(i).spectra(:,:),2);
        [mindb(i).db,idx]=min(s(i).spectra(:,:),[],2);
        [r,c]=ind2sub(sz(1:2:3),idx);
        mindb(i).latlon=s(i).latlon(r,:);
        mindb(i).slow=s(i).slow;
        mindb(i).freq=s(i).freq(c);
        [maxdb(i).db,idx]=max(s(i).spectra(:,:),[],2);
        [r,c]=ind2sub(sz(1:2:3),idx);
        maxdb(i).latlon=s(i).latlon(r,:);
        maxdb(i).slow=s(i).slow;
        maxdb(i).freq=s(i).freq(c);
    else
        % get values for the entire spectra
        sz=size(s(i).spectra);
        meddb(i).db=median(s(i).spectra(:));
        [mindb(i).db,idx]=min(s(i).spectra(:));
        [r,c,p]=ind2sub(sz,idx);
        mindb(i).latlon=s(i).latlon(r,:);
        mindb(i).slow=s(i).slow(c);
        mindb(i).freq=s(i).freq(p);
        [maxdb(i).db,idx]=max(s(i).spectra(:));
        [r,c,p]=ind2sub(sz,idx);
        maxdb(i).latlon=s(i).latlon(r,:);
        maxdb(i).slow=s(i).slow(c);
        maxdb(i).freq=s(i).freq(p);
    end
    
    % deal with f/s averaging
    if(~s(i).vector(1))
        % frequency average
        mindb(i).freq=s(i).freq.';
        maxdb(i).freq=s(i).freq.';
    end
    if(~s(i).vector(2))
        % slowness average
        mindb(i).slow=s(i).slow.';
        maxdb(i).slow=s(i).slow.';
    end
end

end
