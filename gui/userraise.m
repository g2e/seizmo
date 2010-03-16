function [data,rai,fh]=userraise(data,varargin)
%USERRAISE    Interactively scale SEIZMO records by a power
%
%    Usage:    data=userraise(data)
%              data=userraise(data,'field1',value1,...,'fieldN',valueN)
%              [data,ray]=userraise(...)
%              [data,ray,fh]=userraise(...)
%
%    Description: DATA=USERRAISE(DATA) presents an interactive menu and
%     record section plot (arranged by degree distance) to facilitate
%     scaling records in SEIZMO struct DATA by a specified power.  This is
%     useful for noise suppression/enhancement.  The default power is 1
%     which leaves the data as it originally was.  Using power>1 will
%     emphasize high data values while 0<power<1 will emphasize small data
%     values (in an absolute value sense).  Using a power of 0 is
%     equivalent to a 1-bit filter.  Using a power<0 is a bad idea!
%
%     DATA=USERRAISE(DATA,'FIELD1',VALUE1,...,'FIELDN',VALUEN) passes
%     field/value pairs to RECORDSECTION to allow plot customization.
%
%     [DATA,RAI]=USERRAISE(...) returns a struct RAI with the following
%     fields:
%      RAI.power  --  power applied to data
%
%     [DATA,RAI,FH]=USERRAISE(...) returns the record section's figure
%      handle in FH.
%
%    Notes:
%
%    Header changes: DEPMIN, DEPMAX, DEPMEN
%
%    Examples:
%     Scaling records in your dataset may be useful for noise suppression
%     before some other tasks:
%      data=userraise(data);
%      xc=correlate(data,'npeaks',3,'spacing',10);
%      [arr,err,pol]=ttsolve(xc);
%
%    See also: USERWINDOW, USERTAPER, USERMOVEOUT, USERALIGN

%     Version History:
%        Mar. 16, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 16, 2010 at 15:25 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
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

% attempt scaling
try
    % default power
    rai.power=1;
    
    % outer loop - only breaks free on user command
    happy_user=false;
    while(~happy_user)
        % plot records
        fh=recordsection(raise(data,rai.power),varargin{:});
        
        % get choice from user
        pstr=num2str(rai.power);
        choice=menu(...
            ['Alter Power-law Scaling of Data (Currently Power=' ...
            pstr ')?'],'YES','NO');
        
        % act on choice
        if(choice==1)
            % ask user for a particular value
            choice=menu('Scale data by what power?',...
                ['CURRENT (' pstr ')'],...
                '4','3','2','1/2','1/3','1/4','CUSTOM');
            
            switch choice
                case 1 % current
                    % do nothing
                case 2
                    rai.power=4;
                case 3
                    rai.power=3;
                case 4
                    rai.power=2;
                case 5
                    rai.power=1/2;
                case 6
                    rai.power=1/3;
                case 7
                    rai.power=1/4;
                case 8 % custom
                    tmp=inputdlg(...
                        ['Scale data by what power? [' pstr ']:'],...
                        'Custom Power Scaling',1,{pstr});
                    if(~isempty(tmp))
                        try
                            tmp=str2double(tmp{:});
                            if(isscalar(tmp) && isreal(tmp))
                                rai.power=tmp;
                            end
                        catch
                            % do not change rai.power
                        end
                    end
            end
            
            % close old figure
            close(fh(ishandle(fh)));
        else
            happy_user=true;
        end
    end
    
    % apply power
    if(rai.power~=1)
        data=raise(data,rai.power);
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
