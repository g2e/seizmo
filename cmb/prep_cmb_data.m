function [cmt]=prep_cmb_data(indir,outdir,sodcsv,src)
%PREP_CMB_DATA    Prepare core-diffracted dataset
%
%    Usage:    cmt=prep_cmb_data(indir,outdir,sodcsvfile)
%              cmt=prep_cmb_data(indir,outdir,sodcsvfile,datasrc)
%
%    Description:
%     CMT=PREP_CMB_DATA(INDIR,OUTDIR,SODCSVFILE) prepares the SEIZMO data
%     files in INDIR (which contains directories of data corresponding to
%     individual earthquakes).  Preparation includes:
%      ( 1) EVENT SELECTION (by user)
%      ( 2) CMT SELECTION (based on sodcsvfile)
%      ( 3) HEADER CORRECTIONS
%      ( 4) MERGING
%      ( 5) WINNOWING (>600s records only)
%      ( 6) DETRENDING
%      ( 7) TAPERING
%      ( 8) RESAMPLING (1Hz)
%      ( 9) RESPONSE REMOVAL
%      (10) FORCE VERTICALS TO BE ORIENTED UPWARD
%      (11) ROTATION (into radial & transverse)
%      (12) ALIGNMENT (on Pdiff or Sdiff)
%      (13) WINNOWING (require 300s about arrival)
%      (14) RENAMING (rdseed based naming)
%     SODCSVFILE should match the directory list in INDIR.  If it does not
%     then things will fail.  So basically you must use SOD if you want to
%     proceed.  Sorry.  CMT is the GlobalCMT moment tensors used.
%
%     CMT=PREP_CMB_DATA(INDIR,OUTDIR,SODCSVFILE,DATASRC) indicates where
%     the SAC files came from: 'RDSEED' for SAC files output from a SEED
%     volume by rdseed or 'SOD' for SAC files downloaded directly via sod.
%     This is just to call the appropriate header fixing routine.  The
%     default is 'RDSEED'.
%
%    Notes:
%     - If response removal fails, a warning is issued (rather than an
%       error) and the offending dataset is immediately returned rather
%       than the cmts.
%     - INDIR example layout:
%                INDIR
%                 |
%                 +-> EVENTDIR1
%                 +-> EVENTDIR2
%                      |
%                      +-> FILE1
%                      +-> FILE2
%
%    Examples:
%     % You can skip giving a CSV input if you would rather select it
%     % using a GUI menu:
%     cmt=prep_cmb_data('myrawdata','myprepdata');
%
%    See also: CMB_1ST_PASS, CMB_OUTLIERS, CMB_2ND_PASS, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, MAP_CMB_PROFILES

%     Version History:
%        Dec. 12, 2010 - added docs
%        Jan. 12, 2011 - use point_vecticals_upward
%        Jan. 30, 2011 - data source argument, note about response crashes
%        Mar. 19, 2011 - better catches when no data left, deal with
%                        magnitude types that are not in globalcmt list
%        Apr.  4, 2011 - remove response to displacement
%        Feb.  7, 2012 - merge to meld update
%        Mar.  1, 2012 - fixes for scalar structs
%        Mar. 15, 2012 - fixes for pick functions
%        Jan. 28, 2013 - some records have inf data values so some code
%                        was added to catch & delete those records
%        Feb. 26, 2013 - helpful detection of .csv/data desync
%        Apr.  8, 2013 - 2hr desync required due to begin times of some
%                        records being outside 1hr limit
%        Jan. 27, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% default sodcsv (user gui select)
if(nargin==2); sodcsv=[]; end

% default/check
if(nargin<4 || isempty(src)); src='rdseed'; end
if(~isstring(src) || ~any(strcmpi(src,{'rdseed' 'sod'})))
    error('seizmo:prep_cmb_data:badInput',...
        'DATASRC must be either ''RDSEED'' or ''SOD''!');
end

% directory separator
fs=filesep;

% check indir
if(~isstring(indir))
    error('seizmo:prep_cmb_data:badInput',...
        'INDIR must be a string giving one directory!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~isdir(indir))
    error('seizmo:prep_cmb_data:badInput',...
        'INDIR must be a directory!');
end

% check outdir
if(~isstring(outdir))
    error('seizmo:prep_cmb_data:badInput',...
        'OUTDIR must be a string giving one directory!');
end
if(~isabspath(outdir)); outdir=[pwd fs outdir]; end
if(exist(outdir,'file') && ~isdir(outdir))
    error('seizmo:prep_cmb_data:badInput',...
        'OUTDIR must be a directory!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
datelist=char(strcat({dates.name}.'));

% get event list
event=readsodeventcsv(sodcsv);

% get user selected events
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(dates),...
          'ListSize',[170 300],...
          'ListString',datelist);

% silence annoying checkheader warnings
checkoperr('MULTIPLE_DELTA','IGNORE',...
           'INVALID_IZTYPE','IGNORE',...
           'KM_DEPTH','IGNORE');

% loop over user selected events
for i=s(:)'
    % print date
    disp(dates(i).name);
    
    % read in data
    data=readseizmo([indir filesep dates(i).name]);
    
    % clean out data with inf values
    baddata=false(numel(data),1);
    for b=1:numel(data)
        if(any(isinf(data(b).dep))); baddata(b)=true; end
    end
    data(baddata)=[];
    
    % find event cmt
    ievent=ssidx(event,i);
    if(any(strcmpi(ievent.magnitudeType,{'mw' 'mb' 'ms'})))
        cmt(i)=findcmt('time',ievent.time,...
            'location',[ievent.latitude ievent.longitude],...
            'depth',ievent.depth,...
            'magtype',char(ievent.magnitudeType),...
            'magnitude',ievent.magnitude);
    else
        cmt(i)=findcmt('time',ievent.time,...
            'location',[ievent.latitude ievent.longitude],...
            'depth',ievent.depth);
    end
    
    % fix header info
    data=setevent(data,cmt(i));
    switch lower(src)
        case 'rdseed'
            data=fix_rdseed_v48(data);
        case 'sod'
            data=fix_sod_v222(data);
    end
    
    % rough check for major desync between sodcsv & data directories
    % - error if beginning of any record does
    %   not begin within 2 hours of event time
    b=getheader(data,'b');
    if(any(abs(b)>7200))
        error('seizmo:prep_cmb_data:desync',...
            ['Data Directories & .csv file appear desynced!\n' ...
            'Did you forget to send/download a breqfast/SEED file?\n' ...
            'Did you send/download an extra breqfast/SEED file?\n' ...
            'Event: %s\nDirectory: %s'],datestr(ievent.time),...
            dates(i).name);
    end
    
    % merge data
    data=meld(data);
    
    % remove records less than 10 minutes in length
    [b,e]=getheader(data,'b','e');
    data(e-b<600)=[];
    
    % skip to next if none left
    if(isempty(data)); continue; end
    
    % detrend, taper, resample
    data=removetrend(data);
    data=taper(data);
    data=syncrates(data,1,1e-5);
    
    % remove response
    try
        data=removesacpz(data,'tl',1./[200 100]);
    catch
        msg=lasterror;
        warning(msg.identifier,msg.message);
        cmt=data;
        return;
    end
    
    % work with verticals
    vdata=data(vertcmp(data));
    if(~isempty(vdata))
        % rename vertical records to rdseed style
        vdata=genname(vdata,'rdseed');
        
        % flip polarity of verticals if pointing down
        vdata=point_verticals_upward(vdata);
        
        % add arrivals
        vdata=arrivals2picks(vdata,'P,Pdiff,PKP');
        
        % set Pdiff as reference time
        vdata=timeshift(vdata,-findpicks(vdata,'P,Pdiff',1));
        
        % remove records without 300s before & after Pdiff
        [b,e]=getheader(vdata,'b','e');
        vdata(b>-300 | e<300)=[];
    end
    
    % find/rotate horizontals
    rdata=rotate(data);
    if(~isempty(data))
        % rename rotated records
        rdata=genname(rdata,'rdseed');
        
        % add arrivals
        rdata=arrivals2picks(rdata,'S,Sdiff,SKS,sSKS,pSKS');
        
        % set Sdiff as reference time
        rdata=timeshift(rdata,-findpicks(rdata,'S,Sdiff',1));
        
        % remove records without 300s before & after Sdiff
        [b,e]=getheader(rdata,'b','e');
        rdata(b>-300 | e<300)=[];
    end
    
    % save records
    if(isempty([vdata; rdata])); continue; end
    writeseizmo([vdata; rdata],'path',[outdir filesep dates(i).name]);
end

% compact cmt output
cmt=sscat(cmt);

end
