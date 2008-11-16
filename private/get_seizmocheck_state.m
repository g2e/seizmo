function [state]=get_seizmocheck_state()
global SEIZMO
try
    state=SEIZMO.SEIZMOCHECK.ON;
catch
    state=true;
end
end