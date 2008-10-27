function [state]=get_chkhdr_state()
global SACLAB
try
    state=SACLAB.CHKHDR.SKIP;
catch
    state=false;
end
end