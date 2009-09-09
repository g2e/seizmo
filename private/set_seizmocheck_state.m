function []=set_seizmocheck_state(state)
%SET_SEIZMOCHECK_STATE    Turn SEIZMOCHECK on (TRUE) or off (FALSE)
global SEIZMO
SEIZMO.SEIZMOCHECK.ON=state;
end
