function [state]=get_seizchk_state()
global SACLAB
try
    state=SACLAB.SEIZCHK.ON;
catch
    state=true;
end
end