function [state]=get_chkhdr_state()
global SACLAB
try
    state=SACLAB.CHKHDR.ON;
catch
    state=true;
end
end