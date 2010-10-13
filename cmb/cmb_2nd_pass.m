function [results2]=cmb_2nd_pass(results0,bank,gcrng)
%CMB_2ND_PASS    Narrow-band core-diff alignment + normalization
%
%    Usage:    results=cmb_2nd_pass(results,bank,gcrng)
%
%    Notes:
%     bank=filter_bank([0.0125 0.125],'variable',0.2,0.1);
%     gcrng=[90 155];
%     results2=cmb_2nd_pass(results,bank,gcrng);

% todo:

% check results (needs 1st pass, outlier, corrections)
reqfields={'useralign' 'filter' 'usersnr' 'corrections' 'outliers' ...
    'tt_start' 'phase' 'runname' 'dirname' 'cluster' 'usercluster'};
if(~isstruct(results0) || any(~isfield(results0,reqfields)))
    error('seizmo:cmb_2nd_pass:badInput',...
        ['RESULTS must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check filter bank
if(size(bank,2)~=3 ...
        || any(bank(:)<=0 | isnan(bank(:)) | isinf(bank(:))) ...
        || any(bank(:,1)<=bank(:,2) | bank(:,3)<=bank(:,1)))
    error('seizmo:cmb_2nd_pass:badInput',...
        'BANK must be in the format from FILTER_BANK!');
end

% check gcrng
if(~isreal(gcrng) || numel(gcrng)~=2 || gcrng(1)>gcrng(2) ...
        || any(gcrng)<0)
    error('seizmo:cmb_2nd_pass:badInput',...
        'GCRNG must be [MIN MAX] in degree distance!');
end

% loop over each event
cnt=0;
for i=1:numel(results0)
    % run name
    runname=regexprep(results0(i).runname,'1stPass','2ndPass');
    disp(runname);
    
    % skip if no useralign
    if(isempty(results0(i).useralign)); continue; end
    
    % read in data
    data=readseizmo(strcat(results0(i).dirname,filesep,...
        {results0(i).useralign.data.name}'));
    
    % great circle distance
    gcarc=getheader(data,'gcarc');
    
    % align data using 1stPass results
    arr=results0(i).useralign.solution.arr;
    pol=results0(i).useralign.solution.pol;
    data=multiply(data,pol);
    data=timeshift(data,-getheader(data,'o')-arr);
    
    % loop over good clusters
    for j=find(results0(i).usercluster.good(:)')
        % get cluster info
        sj=num2str(j);
        good=find(results0(i).usercluster.T==j ...
            & gcarc>=gcrng(1) & gcarc<gcrng(2));
        
        % census
        pop=numel(good);
        if(pop<3); continue; end
        
        % extract appropriate corrections
        correct=results0(i).corrections;
        correct=fixcorrstruct(correct,good);
        
        % multiband alignment
        results=multibandalign(data(good),bank,[runname '_cluster_' sj],...
            'absxc',false,'estarr',0,'estpol',1,'wgtpow',2);

        % loop over each band in the result to add more info
        for k=1:numel(results)
            % matlab bug workaround
            % - see cmb_1st_pass
            if(isempty(results(k).usersnr))
                results(k).usersnr=[];
            end

            % add run name, quake name
            results(k).phase=results0(i).phase;
            results(k).runname=runname;
            results(k).dirname=results0(i).dirname;

            % add corrections
            results(k).corrections=correct;
        end

        % save results
        save([runname '_cluster_' sj '_results.mat'],'results');

        % export to command line too
        cnt=cnt+1;
        results2(:,cnt)=results;
    end
end

end


function [s]=fixcorrstruct(s,good)
fields=fieldnames(s);
for i=1:numel(fields)
    if(isstruct(s.(fields{i})))
        s.(fields{i})=fixcorrstruct(s.(fields{i}),good);
    else
        s.(fields{i})=s.(fields{i})(good);
    end
end
end

