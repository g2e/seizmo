function [state]=get_checkheader_state()
%GET_CHECKHEADER_STATE    Returns TRUE if CHECKHEADER is on, FALSE if not
global SEIZMO
try
    state=SEIZMO.CHECKHEADER.ON;
catch
    state=true;
end
end