%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure: m-n planes with varied a b %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function AIF_signaling_parameter_screening_main_rescaled_stable(u)

tic

N = 10000;

IDD_matrix = [];
count = 0;
IDD_matrix_temp = AIF_signaling_RecordScore_RandomSampling_rescaled_stable(N,u);
dcount = size(IDD_matrix_temp,1);
IDD_matrix(count+1:count+dcount,:) = IDD_matrix_temp;
count = count + dcount;

IDD_matrix = sortrows(IDD_matrix,3);
filename = sprintf( ...
    "AIF_signaling_parameter_screening_RandomParameter_rescaled_stable_u=%g.xlsx", u);
writematrix(IDD_matrix,filename);
toc

end