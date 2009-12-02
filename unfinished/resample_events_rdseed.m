function []=resample_events_rdseed(indir,outdir)

% get event dirs in indir
dirs=xdir(indir);
dirs={dirs.name}';
dirs(strcmp(dirs,'.') | strcmp(dirs,'..'))=[];

% loop over dirs
for i=1:numel(dirs)
    disp(dirs{i})
    data=readseizmo([indir filesep dirs{i}]);
    data=syncrates(data,1);
    writeseizmo(data,'path',[outdir filesep dirs{i}]);
end