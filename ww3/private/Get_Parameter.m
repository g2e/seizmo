%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch Parameter from Parameter_Table             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameter=Get_Parameter(ival,idx)
global Parameter_Table
parameter=char(Parameter_Table{ival+1}(idx));

