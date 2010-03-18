function []=userclusteredalign(data,varargin)
%USERCLUSTEREDALIGN    Interactive alignment for clustered SEIZMO records
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar. 16, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 16, 2010 at 15:25 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data (dep)
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
end

% attempt align
try
    % basic checks on optional inputs
    if(mod(nargin-1,2))
        error('seizmo:useralign:badInput',...
            'Unpaired OPTION/VALUE!');
    elseif(~iscellstr(varargin(1:2:end)))
        error('seizmo:useralign:badInput',...
            'All OPTIONs must be specified with a string!');
    end
    
    % default correlate options
    info.correlate.npeaks=3;
    info.correlate.spacing=10;
    info.correlate.absxc=true;
    
    % parse correlate options (remove them too)
    keep=true(nargin-1,1);
    for i=1:2:nargin-1
        switch lower(varargin{i})
            case 'npeaks'
                if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                        || varargin{i+1}~=fix(varargin{i+1}) ...
                        || varargin{i+1}<1)
                    error('seizmo:useralign:badInput',...
                        'NPEAKS must be an integer >=1 !');
                end
                info.correlate.npeaks=varargin{i+1};
                keep(i:i+1)=false;
            case 'spacing'
                if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                        || varargin{i+1}<0)
                    error('seizmo:useralign:badInput',...
                        'SPACING must be a positive real (in seconds)!');
                end
                info.correlate.spacing=varargin{i+1};
                keep(i:i+1)=false;
            case 'absxc'
                if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                    error('seizmo:useralign:badInput',...
                        'ABSXC must be a logical value!');
                end
                info.correlate.absxc=varargin{i+1};
                keep(i:i+1)=false;
        end
    end
    varargin=varargin(keep);
    
    % outer loop - only breaks free on user command
    happy_user=false;
    while(~happy_user)
        % usermoveout
        [data0,info.usermoveout,info.figurehandles(1)]=usermoveout(data);

        % userwindow
        [data0,info.userwindow,info.figurehandles(2:3)]=userwindow(data0);

        % usertaper
        [data0,info.usertaper,info.figurehandles(4:5)]=usertaper(data0);
        
        % userraise
        [data0,info.userraise,info.figurehandles(6)]=userraise(data0);

        % menu for correlate options
        while(1)
            % present current settings
            tmp1=num2str(info.correlate.npeaks);
            tmp2=num2str(info.correlate.spacing);
            tmp3='YES'; if(info.correlate.absxc); tmp3='NO'; end
            choice=menu('CHANGE CORRELATE SETTINGS?',...
                ['NUMBER OF PEAKS (' tmp1 ')'],...
                ['PEAK SPACING (' tmp2 's)'],...
                ['ALL POLARITIES ARE MATCHED (' tmp3 ')'],...
                'NO, GO AHEAD AND CORRELATE DATA');

            % proceed by user choice
            switch choice
                case 1 % npeaks
                    choice=menu('NUMBER OF PEAKS TO PICK',...
                        ['CURRENT (' tmp1 ')'],'1','3','5','CUSTOM');
                    switch choice
                        case 1 % CURRENT
                            % leave alone
                        case 2 % 1
                            info.correlate.npeaks=1;
                        case 3 % 3
                            info.correlate.npeaks=3;
                        case 4 % 5
                            info.correlate.npeaks=5;
                        case 5 % CUSTOM
                            tmp=inputdlg(...
                                ['Number of Peaks to Pick? [' tmp1 ']:'],...
                                'Custom Number of Peaks',1,{tmp1});
                            if(~isempty(tmp))
                                try
                                    tmp=str2double(tmp{:});
                                    if(isscalar(tmp) && isreal(tmp) ...
                                            && tmp==fix(tmp) && tmp>=1)
                                        info.correlate.npeaks=tmp;
                                    end
                                catch
                                    % do not change info.correlate.npeaks
                                end
                            end
                    end
                case 2 % spacing
                    tmp=inputdlg(...
                        ['Minimum Spacing Between Peaks (in seconds)? [' ...
                        tmp2 ']:'],'Peak Spacing',1,{tmp2});
                    if(~isempty(tmp))
                        try
                            tmp=str2double(tmp{:});
                            if(isscalar(tmp) && isreal(tmp) && tmp>=0)
                                info.correlate.spacing=tmp;
                            end
                        catch
                            % do not change info.correlate.spacing
                        end
                    end
                case 3 % absxc
                    choice=menu('DO THE POLARITIES ALL MATCH?',...
                        ['CURRENT (' tmp3 ')'],'YES','NO');
                    switch choice
                        case 1 % CURRENT
                            % leave alone
                        case 2 % looking at just positive peaks
                            info.correlate.absxc=false;
                        case 3 % looking at both pos/neg peaks
                            info.correlate.absxc=true;
                    end
                case 4 % all good
                    break;
            end
        end

        % correlate (menu to edit options)
        xc=correlate(data0,...
            'npeaks',info.correlate.npeaks,...
            'spacing',info.correlate.spacing,...
            'absxc',info.correlate.absxc);
        
        % usercluster
        
        % population cut
        
        % fine tuned selection
        % - plotclusters!!!
        % - selectclusters!!!
        
        % loop over groups

        % solve alignment
        [arr,err,pol,zmean,zstd,nc,info.ttsolve]=ttsolve(xc,varargin{:});

        % plot alignment
        info.figurehandles(7)=recordsection(...
            multiply(timeshift(data0,arr),pol));
        
        % arr cut?
        
        % err cut?
        
        % plot alignment again
        
        % end loop over groups
        % - how are we going to export this data
        %   - grp(i).group
        %           .members
        %           .arrival
        %           .stderr
        %           .polarity
        %           .zmean
        %           .zstd
        %           .nc (diagnostic - not useful here)

        % ask user if they are happy with alignment
        choice=menu('KEEP THESE ALIGNMENT(S)?',...
            'YES','NO - TRY AGAIN','NO - CRASH!');
        switch choice
            case 1 % rainbow's end
                happy_user=true;
            case 2 % never never quit!
                close(info.figurehandles(ishandle(info.figurehandles)));
            case 3 % i bear too great a shame to go on
                error('seizmo:useralign:killYourSelf',...
                    'User demanded Seppuku!');
        end
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
