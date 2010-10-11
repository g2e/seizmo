function [results]=cmb_1st_pass(phase,indir)
%CMB_1ST_PASS    Initial core-diff alignment + normalization + clustering
%
%    Usage:    results=cmb_1st_pass(phase,indir)

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check phase
valid={'Pdiff' 'SHdiff' 'SVdiff'};
if(~isstring(phase) || ~ismember(phase,valid))
    error('seizmo:cmb_1st_pass:badPhase',...
        ['PHASE must be one of the following:\n' ...
        sprintf('''%s'' ',valid{:}) '!']);
end

% check indir
if(~isstring(indir))
    error('INDIR must be a string giving one directory!');
elseif(~isdir(indir))
    error('INDIR must be a directory!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
datelist=char(strcat({dates.name}.'));

% get user selected start date
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(dates),...
          'ListSize',[170 300],...
          'ListString',datelist);

% loop over events
for i=1:numel(s)
    % echo directory name to screen
    disp(dates(s(i)).name);
    
    % read in headers
    data=readheader([indir filesep dates(s(i)).name]);
    runname=[dates(s(i)).name '_' phase '_1stPass'];
    
    % appropriate components for this phase
    switch phase
        case 'Pdiff'
            % only vertical data
            data=data(vertcmp(data));
        case 'SHdiff'
            % only transverse component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHT') ...
                | strcmp(kcmpnm,'HHT') ...
                | strcmp(kcmpnm,'LHT'));
        case 'SVdiff'
            % only radial component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHR') ...
                | strcmp(kcmpnm,'HHR') ...
                | strcmp(kcmpnm,'LHR'));
    end
    
    % read data
    data=readdata(data);
    
    % sort by distance
    data=sortbyfield(data,'gcarc');
    
    % set up a singleband for multibandalign
    bank=1./[42.5 80 25];
    
    % multibandalign
    tmp=multibandalign(data,bank,runname,'absxc',true,'estarr',[],'wgtpow',2);
    
    % matlab bug workaround
    % - really strange...saving fails to write all fields if we dont do this
    if(isempty(tmp.usersnr))
        tmp.usersnr=[];
    end
    
    % add run name, quake name
    tmp.phase=phase;
    tmp.runname=runname;
    tmp.dirname=[indir filesep dates(s(i)).name];
    
    % add corrections
    if(~isempty(tmp.useralign))
        tmp.corrections=cmb_corrections(phase,tmp.useralign.data);
    else
        tmp.corrections=[];
    end
    
    % save results
    save([runname '_results.mat'],'-struct','tmp');
    
    % read in to check
    tmp0=load([runname '_results.mat']);
    reqfields={'useralign' 'filter' 'usersnr' 'tt_start' 'phase' ...
        'runname' 'dirname' 'corrections'};
    if(any(~isfield(tmp0,reqfields)))
        disp(runname);
        error('seizmo:cmb_1st_pass:badWrite',...
            'Saving results to file failed!')
    end
    
    % export to command line too
    results(i)=tmp;
end

end
