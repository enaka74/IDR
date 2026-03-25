%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure: m-n planes with varied a b %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constants
param_range_1 = [1,2];
param_range_2 = [1,2];
%param_range_1 = logspace(1,1,1);
%param_range_2 = logspace(2,2,1);

f = figure(1);
f.Position= [0, 0, 500, 500];
tiledlayout(length(param_range_2),length(param_range_1),"TileSpacing","tight");

[param_list_1,param_list_2] = meshgrid(param_range_1,param_range_2);
parfor i=1:numel(param_list_1)
    IDD_create_mnMap_AUC_axay(i,[param_list_1(i),param_list_2(i)])
end

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure: m-n planes with varied a b %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Constants
param_range_1 = logspace(-1,2,4);
param_range_2 = logspace(-1,2,4);

f = figure(1);
f.Position= [0, 0, 600, 1000];
tiledlayout(length(param_range_2),length(param_range_1),"TileSpacing","loose");

[param_list_1,param_list_2] = meshgrid(param_range_1,param_range_2);
parfor i=1:numel(param_list_1)
    IDD_create_mnMap(i,[param_list_1(i),param_list_2(i)])
end
%}