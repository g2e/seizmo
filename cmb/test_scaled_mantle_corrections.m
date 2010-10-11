function [varred]=test_scaled_mantle_corrections(results,scale)
%TEST_SCALED_MANTLE_CORRECTIONS    Tests variance reduction of scaled mantle corrections
%
%    Usage:    varred=test_scaled_mantle_corrections(results,scales)

% - we want to find the best multipliers for the corrections given here
%   - multiplier that best eliminates residuals from linear fit
%   - each cluster separately
%   - weight by cluster size (n-2)
%   - variance reduction = (1-(after/before))*100

% check nargin

% check results (needs 1st pass, outlier, corrections)
reqfields={'useralign' 'filter' 'usersnr' 'corrections' 'outliers' ...
    'tt_start' 'phase' 'runname' 'dirname' 'cluster' 'usercluster'};
if(~isstruct(results) || any(~isfield(results,reqfields)))
    error('seizmo:cmb_outliers:badInput',...
        ['RESULTS must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check scale
if(~isreal(scale))
    error('seizmo:test_scaled_mantle_corrections:badInput',...
        'SCALE must be a real-valued array!');
end

% preallocate variance reduction matrix
varred=zeros(size(scale));
varc=varred;

% loop over each event
rcnt=0;
for i=1:numel(results)
    % loop over each good cluster
    for j=find(results(i).usercluster.good(:)')
        % get cluster info
        good=find(results(i).usercluster.T==j & ~results(i).outliers);
        pop=numel(good);
        
        % just to make sure
        if(pop<=2); continue; end
        
        % get variance about fit before correction
        dd=getheader(results(i).useralign.data(good),'gcarc');
        arr=results(i).useralign.solution.arr(good);
        
        % get variance about linear fit before correction
        m=wlinem(dd,arr)';
        resid=arr-polyval(fliplr(m),dd);
        var0=var(resid);
        
        % get variance reduction after corrections
        for k=1:numel(scale)
            arrcorr=arr-(results(i).corrections.ellcor(good)...
                +results(i).corrections.crucor.ak135(good)...
                +scale(k)*results(i).corrections.mancor.mitp08.full(good));
            m=wlinem(dd,arrcorr)';
            resid=arrcorr-polyval(fliplr(m),dd);
            varc(k)=var(resid);
            if(isnan(var(resid)))
                disp([i j]);
            end
            varred(k)=varred(k)+(pop-2)*(1-(var(resid)/var0))*100;
        end
        
        % increment counters
        rcnt=rcnt+pop-2;
    end
end

% divide by weighting to get true variance reduction
varred=varred/rcnt;

end
