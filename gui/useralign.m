function [db,data]=useralign(data,phase,db)
%USERALIGN    Interactively align SEIZMO records
%
%    Usage:    db=useralign(data,phase)
%              db=useralign(data,phase,db)
%
%    Description: DB=USERALIGN(DATA,PHASE) presents a set of interactive
%     menus and plots to facilitate aligning a seismic phase with a few
%     mouse clicks.  The output is a struct containing all relevant
%     parameters utilized in the alignment as well as all alignment info
%     produced.  The layout of struct DB is as follows:
%
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
%    See also: userwindow, usertaper, usercluster

%     Version History:
%        Sep. 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2009 at 07:30 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
disp('CHECKING DATASET')
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% get header setup
[h,idx]=versioninfo(data);
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% number of records
nrecs=numel(data);

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

% require event info is equal
if(size(unique(ev,'rows'),1)>1 || numel(unique(o))>1)
    error('seizmo:useralign:badHeader','Conflicting Event Information!');
end

% synchronize/clean up (origin time precision always >0.001)
disp('SYNCHRONIZING RECORDS TO ORIGIN')
data=timeshift(data,-o,'io');
data=changeheader(data,'o',0);

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


% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
