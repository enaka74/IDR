%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure: m-n planes with varied a b %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function NFL_signaling_parameter_screening_main_minscale()

tic
N = 5000;
IDD_matrix = [];
count = 0;
IDD_matrix_temp = NFL_signaling_RecordScore_RandomSampling_minscale(N);
dcount = size(IDD_matrix_temp,1);
IDD_matrix(count+1:count+dcount,:) = IDD_matrix_temp;
count = count + dcount;

IDD_matrix = sortrows(IDD_matrix,1);
filename = "NFL_signaling_parameter_screening_RandomParameter_minscale.xlsx";
writematrix(IDD_matrix,filename);

toc

end