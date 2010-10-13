function [vol]=geofkarfcorr(vol,mfarfflag)
%GEOFKARFCORR    Computes correlation with ARFs in the parameter space
%
%    Usage:    vol=geofkarfcorr(vol)
%              vol=geofkarfcorr(vol,mfarfflag)
%
%    Description:
%     VOL=GEOFKARFCORR(VOL) correlates the geofk volume VOL with ARFs
%     generated at each point, slowness and frequency in the parameter
%     space of VOL.  The correlation output is placed in the .beam field of
%     VOL and .normdb is set to 0.  This allows using 
%
%     VOL=GEOFKARFCORR(VOL,MFARFFLAG) is for multi-frequency ARFs in the
%     correlation.  This only applies if VOL is averaged across frequency.
%     The default is FALSE and does not use multi-frequency ARFs.
%
%    Notes:
%     - This takes a long time!
%
%    Examples:
%     %
%
%    See also: GEOFKCORR, GEOFKARF, PLOTGEOFKMAP

%     Version History:
%        Oct. 11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default mfarfflag
if(nargin<2 || isempty(mfarfflag)); mfarfflag=false; end

% check input
error(chkgeofkstruct(vol));
if(~islogical(mfarfflag) || ~isscalar(mfarfflag))
    error('seizmo:geofkarfcorr:badInput',...
        'MFARFFLAG must be TRUE or FALSE!');
end

% add mfarfflag
[vol.mfarfflag]=deal(mfarfflag);

% loop over each volume
for i=1:numel(vol)
    % parameter space
    method=vol(i).method;
    stlalo(:,1)=vol(i).stla;
    stlalo(:,2)=vol(i).stlo;
    latlon=vol(i).latlon;
    horzslow=vol(i).horzslow;
    freq=vol(i).freq;
    
    % check if averaged across frequency or slowness
    % vflag(1) == horzslow
    %      (2) == freq
    vflag=vol(i).volume;
    
    % loop over parameters, generating ARFs and correlating with data
    nll=size(latlon,1);
    ns=numel(horzslow);
    nf=numel(freq);
    if(~vflag(2) && mfarfflag)
        r=nan(nll,ns,1,'single');
    else
        r=nan(nll,ns,nf,'single');
    end
    for l=1:size(latlon,1)
        for s=1:numel(horzslow)
            % frequency averaged?
            if(vflag(2)) % not averaged in frequency
                for f=1:numel(freq)
                    % generate ARF
                    arf=geofkarf(stlalo,latlon,horzslow,latlon(l,:),...
                        horzslow(s),freq(f),method);
                    
                    % horzslow average?
                    if(~vflag(1))
                        arf=geofkarf2map(arf);
                    end
                    
                    % correlate
                    r(l,s,f)=geofkcorr(...
                        geofksubvol(vol,[freq(f) freq(f)]),arf);
                end
            else % averaged in frequency
                % allow multi-freq ARF
                if(mfarfflag)
                    % multi-freq ARF at same point & slowness
                    arf=geofkarf(stlalo,latlon,horzslow,latlon(l,:),...
                        horzslow(s),freq,method);
                    
                    % horzslow average?
                    if(~vflag(1))
                        arf=geofkarf2map(arf);
                    end
                    
                    % correlate
                    r(l,s,1)=geofkcorr(vol,arf);
                else
                    % loop over each frequency
                    for f=1:numel(freq)
                        % generate ARF
                        arf=geofkarf(stlalo,latlon,horzslow,latlon(l,:),...
                            horzslow(s),freq(f),method);
                        
                        % horzslow average?
                        if(~vflag(1))
                            arf=geofkarf2map(arf);
                        end
                        
                        % correlate
                        r(l,s,f)=geofkcorr(vol,arf);
                    end
                end
            end
        end
    end
    
    % replace beam with r
    % - volume needs changing ([1 0/1])
    % - normdb == 0
    vol(i).normdb=0;
    vol(i).volume=[true (~mfarfflag | vflag(2))];
    vol(i).beam=r;
end

end
