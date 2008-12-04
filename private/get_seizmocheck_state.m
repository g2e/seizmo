function [state]=get_seizmocheck_state()
%GET_SEIZMOCHECK_STATE    Returns TRUE if SEIZMOCHECK is on, FALSE if not
global SEIZMO
try
    state=SEIZMO.SEIZMOCHECK.ON;
catch
    state=true;
end
end