function [state]=get_checkheader_state()
global SEIZMO
try
    state=SEIZMO.CHECKHEADER.ON;
catch
    state=true;
end
end