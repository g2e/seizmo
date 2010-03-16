function [db,data]=useralign(data,varargin)
%USERALIGN    Interactively align SEIZMO records
%
%    Usage:    db=useralign(data)
%              db=useralign(data,'option1',value,...,'optionN',value)
%
%    Description: DB=USERALIGN(DATA,PHASE) presents a set of interactive
%     menus and plots to facilitate aligning a seismic phase with a few
%     mouse clicks.  The output is a struct containing all relevant
%     parameters utilized in the alignment as well as all alignment info
%     produced.  The layout of struct DB is as follows:
%      DB.run.time
%            .event_idx
%            .phase_idx
%            .bank_idx
%            .filter_idx
%            .taper.type
%                  .width
%                  .option
%            .snr.cutoff
%                .cut
%            .correct.for_mantle
%                    .for_crust
%                    .for_elev
%                    .for_ellip
%                    .for_previous_run
%            .arrinv.npeaks
%                   .abspeaks
%                   .maxiter
%                   .results
%            .cluster.dissim_cutoff
%                    .groups
%                    .pop_cutoff
%            .arrival.station_idx
%                    .cmpinc
%                    .cmpaz
%                    .correct.mantle_idx
%                            .crust_idx
%                            .elev_idx
%                            .ellip_idx
%                            .run_idx
%                    .signal_window
%                    .noise_window
%                    .ph_snr
%                    .ph_arr
%                    .ph_arr_err
%                    .ph_amp
%                    .ph_amp_err
%                    .ph_z_mean
%                    .ph_z_std
%                    .env_snr
%                    .env_arr
%                    .env_arr_err
%                    .env_amp
%                    .env_amp_err
%                    .env_z_mean
%                    .env_z_std
%            .slow.arrival_idx1
%                 .arrival_idx2
%                 .ph_slowness
%                 .ph_slowness_err
%                 .ph_decay
%                 .ph_decay_err
%                 .env_slowness
%                 .env_slowness_err
%                 .env_decay
%                 .env_decay_err
%      DB.event.ev
%              .nz
%      DB.phase.name
%              .model1d_idx
%      DB.bank.type
%             .range
%             .width
%             .offset
%             .filter.type
%                    .style
%                    .corners
%                    .order
%      DB.station.knetwk
%                .kstnm
%                .khole
%                .kcmpnm
%                .st
%      DB.model1d.name
%      DB.model3d.name
%                .model1d_idx
%      DB.correct.mantle.model3d_idx
%                       .station_idx
%                       .event_idx
%                       .phase_idx
%                       .correction
%                .crust.model3d_idx
%                      .event_idx
%                      .station_idx
%                      .phase_idx
%                      .model
%                      .correction
%                .elev.model3d_idx
%                     .event_idx
%                     .station_idx
%                     .phase_idx
%                     .elevation
%                     .correction
%                .ellip.model1d_idx
%                      .event_idx
%                      .station_idx
%                      .phase_idx
%                      .correction
%
%     DB=USERALIGN(DATA,PHASE,DB) adds all parameters and alignment info to
%     an already existing struct DB.  This allows centralizing information
%     for several events/bandpasses/runs in one location.
%
%    Notes:
%
%    Examples:
%     Align a suite of bandpasses:
%      db=[];
%      bank=filter_bank([0.01 0.2],'variable',0.4,0.2)
%      for i=1:size(bank,1)
%          data=iirfilter(data,'bandpass','butter',bank(i,2:3),4);
%          db=useralign(data,phase,db);
%      end
%
%    See also: USERWINDOW, USERTAPER, USERCLUSTER

%     Version History:
%        Sep. 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2009 at 07:30 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
[h,idx]=versioninfo(data,'dep');

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undefstr=getsubfield(h,'undef','stype').';
undef=undef(idx);
undefstr=undefstr(idx);

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

% get global
global SEIZMO

% attempt rest
try
    % number of records
    nrecs=numel(data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% HANDLING OPTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % valid values for options
    valid.OPTIONS={'model1d' 'model3d' 'phase' 'correct'};
    valid.MODEL1D={'PREM' 'AK135' 'IASP91'};
    valid.MODEL3D={'HMSL-06' 'PRI-05'};
    valid.PHASE={'P' 'Pdiff' 'S' 'Sdiff'};
    
    % option defaults
    option.DB=[];      % db
    option.MODEL1D=[]; % prem/ak135/iasp91
    option.MODEL3D=[]; % hmsl-06/pri-05
    option.PHASE=[];   % P/Pdiff/S/Sdiff/etc...
    option.CORRECT=[]; % true/false
    
    % id function
    me=mfilename;
    pkgme=['seizmo:' me ':'];
    
    % get options from SEIZMO global
    ME=upper(me);
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
            error([pkgme 'badInput'],...
                'Options must be specified as a string!');
        end
        if(~isempty(varargin{i+1}))
            option.(upper(varargin{i}))=varargin{i+1};
        end
    end
    
    % check options
    fields=fieldnames(option);
    for i=1:numel(fields)
        % get value of field
        value=option.(fields{i});
        
        % check for unknown options
        if(~any(strcmpi(fields{i},valid.OPTIONS)))
            error([pkgme 'badInput'],'Unknown option: %s !',fields{i});
        end
        
        % skip if empty
        if(isempty(value)); continue; end
        
        % specific checks
        switch lower(fields{i})
            case {'model1d' 'model3d'}
                if(~ischar(value) || size(value,1)~=1 ...
                        || ~any(strcmpi(value,valid.(fields{i}))))
                    error([pkgme 'badInput'],...
                        ['%s option must be one of the following:\n'...
                        sprintf('%s ',valid.(fields{i}){:})],fields{i});
                end
            case 'phase'
                if(~ischar(value) || size(value,1)~=1)
                    error([pkgme 'badInput'],...
                        'PHASE must be a string!');
                end
            case 'correct'
                if(~islogical(value) || ~isscalar(value))
                    error([pkgme 'badInput'],...
                        '%s option must be a scalar logical!',fields{i});
                end
            case 'db'
                if(~isstruct(value) || ~isscalar(value))
                    error([pkgme 'badInput'],...
                        '%s option must be a scalar struct!',fields{i});
                end
            otherwise
                % we should check all options
                error([pkgme 'badDeveloper'],...
                    'Failed to check option: %s !',fields{i});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% BUILDING DATABASE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize db
    if(isempty(option.DB))
        db.event([],1)=struct('ev',[],'nz',[]);
        db.station([],1)=struct('knetwk',[],'kstnm',[],'khole',[],...
            'kcmpnm',[],'st',[]);
        db.model1d([],1)=struct('name',[]);
        db.phase([],1)=struct('name',[],'model1d_idx',[]);
        db.model3d([],1)=struct('name',[],'model1d_idx',[]);
        db.bank([],1)=struct('type',[],'range',[],'width',[],...
            'offset',[],'filter',[]);
        db.correct([],1)=struct('raypath',[],'mantle',[],'crust',[],...
            'elev',[],'ellip',[]);
    else
        db=option.DB;
    end
    
    % deal with event info
    [ev,o]=getheader(data,'ev','o');
    if(any(any(ev==undef(:,[1 1 1 1]) | isnan(ev) | isinf(ev))))
        error([pkgme 'badHeader'],...
            'EVLA EVLO EVEL EVDP fields must be defined!');
    end
    if(any(any(o==undef | isnan(o) | isinf(o))))
        error([pkgme 'badHeader'],...
            'O field must be defined!');
    end
    data=timeshift(data,-o,'io');
    nz=getheader(data,'nz');
    if(any(any(nz==undef(:,[1 1 1 1 1 1]) | isnan(nz) | isinf(nz))))
        error([pkgme 'badHeader'],...
            'NZ fields must be defined!');
    end
    if(size(unique(ev,'rows'),1)>1 || size(unique(nz,'rows'),1)>1)
        error([pkgme 'badHeader'],'Conflicting Event Information!');
    end
    [db,event_idx]=db_update(db,'event',struct('ev',ev(1,:),'nz',nz(1,:)));

    % deal with station info
    [st,name]=getheader(data,'st','kname');
    if(any(any(st==undef(:,[1 1 1 1]) | isnan(st) | isinf(st))))
        error([pkgme 'badHeader'],...
            'STLA STLO STEL STDP must be defined!');
    end
    if(any(strcmpi(name,undefstr(:,[1 1 1 1]))))
        error([pkgme 'badHeader'],...
            'KNETWK KSTNM KHOLE KCMPNM must be defined!');
    end
    if(size(unique(name,'rows'),1)~=nrecs)
        error([pkgme 'badHeader'],...
            ['KNETWK KSTNM KHOLE KCMPNM must form a unique\n' ...
            'combination for each record!']);
    end
    station_idx=nan(nrecs,1);
    for i=1:nrecs
        [db,station_idx(i)]=db_update(db,'station',struct('st',st(i,:),...
            'knetwk',name(i,1),'kstnm',name(i,2),'khole',name(i,3),...
            'kcmpnm',name(i,4)));
    end
    
    % ask for phase if none given
    while(isempty(option.PHASE))
        option.PHASE=listdlg('PromptString','Select a phase:',...
            'Name','Seismic Phase',...
            'SelectionMode','single',...
            'ListString',valid.PHASE);
    end
    
    % ask for model1d if none given
    while(isempty(option.MODEL1D))
        option.MODEL1D=listdlg('PromptString','Select a 1d model:',...
            'Name','Seismic 1D Model',...
            'SelectionMode','single',...
            'ListString',valid.MODEL1D);
    end
    [db,model1d_idx]=db_update(db,'model1d',struct('name',option.MODEL1D));
    [db,phase_idx]=db_update(db,'phase',struct('name',option.PHASE,...
        'model1d_idx',model1d_idx));
    
    % ask to correct for 3D effects
    % - maybe this should allow for more fine-grained control
    %   (ie allow only certain corrections)
    while(isempty(option.CORRECT))
        choice=menu('Correct for 3D Structure?','YES','NO');
        if(choice==1)
            option.CORRECT=true;
        elseif(choice==2)
            option.CORRECT=false;
        end
    end
    db.run.correct.for_mantle=option.CORRECT;
    db.run.correct.for_crust=option.CORRECT;
    db.run.correct.for_elev=option.CORRECT;
    db.run.correct.for_ellip=option.CORRECT;
    db.run.correct.for_geospread=option.CORRECT;
    
    if(option.CORRECT)
        % currect hack
        % - loads 3d corrections from data headers
        %   - resp0 = mantle correction
        %   - resp1 = crust correction
        %   - resp2 = elev correction
        %   - resp3 = ellip correction
        %   - resp4 = custom correction
        %   - resp5 = saved correction from previous run
        %
        % - models are 'unknown'
        % - other idx are straight forward
        [mancor,crucor,elevcor,ellipcor,customcorr,prevcorr]=...
            getheader(data,'resp0','resp1','resp2','resp3','resp4','resp5');
        
        % ask for model3d if none given
        %while(isempty(option.MODEL3D))
        %    option.MODEL1D=listdlg('PromptString','Select a 3d model:',...
        %        'Name','Seismic 3D Model',...
        %        'SelectionMode','single',...
        %        'ListString',valid.MODEL3D);
        %end
        %db.model3d.name=option.MODEL3D;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ask user what they want to do
    analysis=menu('Which analysis type?','Single Band','Multi-Band');
    switch analysis
        case 1 % single band
            % save corrections
            % save group info
        case 2 % multi-band
            % no clustering
    end
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end




% 

% get header info
% iftype, leven, ncmp
[st,ev,o,gcarc]=getheader(data,'st','ev','o','gcarc');

% check if undefined
if(any(any(st==undef(:,[1 1 1 1]) | isnan(st) | isinf(st))))
    error('seizmo:useralign:badHeader',...
        'STLA STLO STEL STDP must be defined!');
elseif(any(any(ev(:,[1 2 4])==undef(:,[1 1 1]) ...
        | isnan(ev(:,[1 2 4])) | isinf(ev(:,[1 2 4])))))
    error('seizmo:useralign:badHeader',...
        'EVLA EVLO EVEL EVDP must be defined!');
elseif(any(o==undef | isnan(o) | isinf(o)))
    error('seizmo:useralign:badHeader',...
        'O field must be defined!');
elseif(any(gcarc==undef | isnan(gcarc) | isinf(gcarc)))
    error('seizmo:useralign:badHeader',...
        'GCARC field must be defined!');
end




% get 1d model from user
models_1d={'PREM' 'AK135' 'IASP91'};
m1d=lower(models_1d{menu('CHOOSE A 1D MODEL',models_1d{:})});

% get phase times/paths
% *** this is a quick hack for what addarrivals will do ***
disp('GETTING PREDICTED PHASE ARRIVAL TIMES')
print_time_left(0,nrecs,true);
for i=1:nrecs
    data(i).misc.arrivals=tauppath('m',m1d,'p',phase,'d',gcarc(i),...
        'z',ev(i,4)/1000);
    if(isempty(data(i).misc.arrivals))
        error('seizmo:useralign:noPhaseHere',...
            'PHASE %s does not exist at station/record %d!',phase,i);
    end
    print_time_left(i,nrecs);
end

% create/check db
if(nargin==2 || isempty(db))
    db.event=struct('lat',ev(1,1),'lon',ev(1,2),...
        'elev',ev(1,3),'depth',ev(1,4));
    db.station=struct('lat',[],'lon',[],'elev',[],'depth',[]);
    db.record=struct('path',[],'name',[]);
    db.run=struct('time',datestr(now));
    db.phase=struct('name',[],'model_1d',m1d);
    db.snrcut=struct('noise',[],'signal',[],'cutoff',[]);
    db.correct=struct('event',[],'phase',[],'station',[],...
        'model_3dp',[],'model_3ds',[],'mantlecorr',[]);
    db.cluster=struct('group',[]);
    db.align=struct('event',1);
else
    % check db
    
end

% get travel time corrections
% - have two options
%   - quick way: lookup value in table (computed elsewhere/previously)
%   - slow way: trace path through a velocity model
% *** no tracing yet ***
% *** we need tauppath to give p/s for each path segment ***
models_3dp={'HMSL-P06' 'PRI-P05' 'MIT-P08' 'MIT-P06' 'DZ04'};
models_3ds={'HMSL-S06' 'PRI-S05' 'SB4L18' 'S20RTSb' 'SAW24B16' 'SG2006e'};
m3dp=models_3dp{menu('CHOOSE A 3D P VELOCITY MODEL',models_3dp{:})};
m3ds=models_3ds{menu('CHOOSE A 3D S VELOCITY MODEL',models_3ds{:})};
mantlecorr=[];
crustcorr=[];
ellipcorr=[];
elevcorr=[];

% apply corrections
% - ask user if this is desired
if(menu('APPLY 3D CORRECTIONS?','YES','NO')==1)
    tt=tt+mantlecorr+crustcorr+ellipcorr+elevcorr;
end

% now do a quick snr cut (aks user for snr)
% - usersnrcut
% - histogram or scatter

% user window
[data,win]=userwindow(data,true,@removemean);

% user taper
[data,tpr]=usertaper(data);

% get correlate parameters
% - ask npeaks, absolute, cube, adjacent
prompt={'NUMBER OF CORRELOGRAM'; 'PEAKS TO CONSIDER?'};
npeaks=2*menu(prompt,'1','3','5','7','9')-1;
absxc=[true false];
absxc=absxc(menu('PICK ABSOLUTE VALUE PEAKS?','YES','NO'));
cube=[true false];
cube=cube(menu('CUBE RECORDS?','YES','NO'));
data=seizmofun(data,@(x)x.^3);
% *** skipping adjacent for now ***

% correlate
peaks=correlate(data,'npeaks',npeaks,'absxc',absxc);

% ask to cluster
if(menu('CLUSTER DATA?','YES','NO')==1)
    grp=usercluster(data,peaks.cg(:,:,1,ceil(end/2)));
else
    grp.cutoff=nan;
    grp.method='nan';
    grp.criterion='nan';
    grp.perm=(1:20).';
    grp.color=[1 1 1];
    grp.T=ones(nrecs,1);
end

% invert


end
