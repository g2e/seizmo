function []=set_seizmocheck_state(state)
%SET_SEIZMOCHECK_STATE    Turn SEIZMOCHECK on (true) or off (false)
global SEIZMO
SEIZMO.SEIZMOCHECK.ON=state;
end