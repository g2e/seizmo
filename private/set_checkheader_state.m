function []=set_checkheader_state(state)
%SET_CHECKHEADER_STATE    Turn CHECKHEADER on (true) or off (false)
global SEIZMO
SEIZMO.CHECKHEADER.ON=state;
end