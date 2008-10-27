function [state]=get_seischk_state()
global SACLAB
try
    state=SACLAB.SEISCHK.SKIP;
catch
    state=false;
end
end