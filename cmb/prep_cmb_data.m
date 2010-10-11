function [cmt]=prep_cmb_data(indir,outdir,sodcsv)
%PREP_CMB_DATA    Prepare Core-Diffracted Dataset
%
%    Usage:    cmt=prep_cmb_data(indir,outdir)
%              cmt=prep_cmb_data(indir,outdir,sodcsvfile)

% check nargin
error(nargchk(2,3,nargin));

% default sodcsv (user gui select)
if(nargin==2); sodcsv=[]; end

% check indir
if(~isstring(indir))
    error('INDIR must be a string giving one directory!');
elseif(~isdir(indir))
    error('INDIR must be a directory!');
end

% check outdir
if(~isstring(outdir))
    error('OUTDIR must be a string giving one directory!');
elseif(exist(outdir,'file') && ~isdir(outdir))
    error('OUTDIR must be a directory!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
datelist=char(strcat({dates.name}.'));

% get event list
event=readsodeventcsv(sodcsv);

% get user selected start date
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(dates),...
          'ListSize',[170 300],...
          'ListString',datelist);

% silence checkheader annoying warnings
checkoperr('MULTIPLE_DELTA','IGNORE',...
           'INVALID_IZTYPE','IGNORE');

% loop over user selected events
for i=s(:)'
    % print date
    disp(dates(i).name);
    
    % read in data
    data=readseizmo([indir filesep dates(i).name]);
    
    % find event csv
    cmt(i)=findcmt('time',event(i).time,...
        'location',[event(i).latitude event(i).longitude],...
        'depth',event(i).depth,...
        'magtype',event(i).magnitudeType,...
        'magnitude',event(i).magnitude);
    
    % fix header info
    data=setevent(data,cmt(i));
    data=fix_rdseed_v48(data);
    
    % merge data
    data=merge(data);
    
    % remove records less than 10 minutes in length
    [b,e]=getheader(data,'b','e');
    data(e-b<600)=[];
    
    % detrend, taper, resample
    data=removetrend(data);
    data=taper(data);
    data=syncrates(data,1,1e-5);
    
    % remove response
    try
        data=removesacpz(data,'f',1./[200 100],'units','vel');
    catch
        msg=lasterror;
        warning(msg.identifier,msg.message);
        cmt=data;
        return;
    end
    
    % work with verticals
    vdata=data(vertcmp(data));
    
    % add arrivals
    vdata=addarrivals(vdata,'ph','P,Pdiff,PKP');
    
    % set Pdiff as reference time
    vdata=timeshift(vdata,-getarrival(vdata,{'P' 'Pdiff'}));
    
    % remove records without 300s before & after Pdiff
    [b,e]=getheader(vdata,'b','e');
    vdata(b>-300 | e<300)=[];
    
    % rotate horizontals
    rdata=rotate(data);
    
    % add arrivals
    rdata=addarrivals(rdata,'ph','S,Sdiff,SKS,sSKS,pSKS');
    
    % set Sdiff as reference time
    rdata=timeshift(rdata,-getarrival(rdata,{'S' 'Sdiff'}));
    
    % remove records without 300s before & after Sdiff
    [b,e]=getheader(rdata,'b','e');
    rdata(b>-300 | e<300)=[];
    
    % rename rotated records
    rdata=genname(rdata,'rdseed');
    
    % save records
    writeseizmo([vdata; rdata],'path',[outdir filesep dates(i).name]);
end

end
