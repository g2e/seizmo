function [data]=prep_events_rdseed(indir,outdir)

% get event dirs in indir
dirs=xdir(indir);
dirs={dirs.name}';
dirs(strcmp(dirs,'.') | strcmp(dirs,'..'))=[];

% loop over dirs
for i=116:116
    disp(dirs{i})
    data=merge(fix_rdseed_v48(readseizmo([indir filesep dirs{i}])));
    data=syncrates(...
        removetrend(...
        data((getheader(data,'e')-getheader(data,'b'))>300)),1);
    writeseizmo(data,'path',[outdir filesep dirs{i}]);
end