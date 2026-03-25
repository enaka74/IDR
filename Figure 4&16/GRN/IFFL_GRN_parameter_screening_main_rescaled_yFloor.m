%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure: m-n planes with varied a b %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function IFFL_GRN_parameter_screening_main_rescaled_yFloor(u)

tic

N = 10000;
%N = 284;

IDD_matrix = [];
count = 0;
IDD_matrix_temp = IFFL_GRN_RecordScore_RandomSampling_rescaled_yFloor(N,u);
dcount = size(IDD_matrix_temp,1);
IDD_matrix(count+1:count+dcount,:) = IDD_matrix_temp;
count = count + dcount;

IDD_matrix = sortrows(IDD_matrix,1);
filename = sprintf( ...
    "IFFL_GRN_parameter_screening_RandomParameter_rescaled_yFloor_u=%g.xlsx", u);
writematrix(IDD_matrix,filename);

toc

end