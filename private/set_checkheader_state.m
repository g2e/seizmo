function []=set_checkheader_state(state)
%SET_CHECKHEADER_STATE    Turn CHECKHEADER on (TRUE) or off (FALSE)
global SEIZMO
SEIZMO.CHECKHEADER.ON=state;
end
