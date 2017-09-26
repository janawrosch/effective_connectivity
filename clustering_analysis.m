function [clustering_parameter]=clustering_analysis(correlation_matrix, cutoffzeile, correlation_matrix_bin, show_figures)

%% 1. Clustering coefficient weighted
cc_values=clustering_coef_wd(correlation_matrix);
cc_values_clean=cc_values(cc_values>0);
if isempty(cc_values_clean)==1
    cc_values_clean=0;
end
cc_mean_weighted=mean(cc_values_clean);
cc_std_weighted=std(cc_values_clean);
[cc_wei_n, cc_wei_x]=hist(cc_values_clean, ceil(cutoffzeile/3));
cc_wei_h=kstest(cc_wei_n);


%% 2. Clustering coefficient binary
cc_values_bin=clustering_coef_bd(correlation_matrix_bin);
cc_values_bin=cc_values_bin(cc_values_bin>0);
if isempty(cc_values_bin)==1
    cc_values_bin=0;
end
cc_bin=mean(cc_values_bin);
cc_bin_std=std(cc_values_bin);
[cc_bin_n, cc_bin_x]=hist(cc_values_bin, ceil(cutoffzeile/3));
cc_bin_h=kstest(cc_bin_n);

%% 3. Transitivity
transitivity=transitivity_wd(correlation_matrix);

%% 4. Local efficiency
local_efficiency=efficiency_wei(correlation_matrix,1);
local_efficiency_clean=local_efficiency(local_efficiency>0);
if isempty(local_efficiency_clean)==1
    local_efficiency_clean=0;
end
local_efficiency_mean=mean(local_efficiency_clean);
local_efficiency_std=std(local_efficiency_clean);

[le_n, le_x]=hist(local_efficiency_clean, ceil(cutoffzeile/3));
le_h=kstest(le_n);

%% 5. Modularity

[~, modularity_degree]=modularity_dir(correlation_matrix);

%% 6. Show figures
if show_figures==1
    figure
    hist(cc_values_clean,10)
    xlabel('Clustering coefficient - weighted')
    ylabel('Count')
    title(sprintf('Clustering coefficient mean: %f; SD: %f', cc_mean_weighted, cc_std_weighted))

    figure
    hist(cc_values_bin,10)
    xlabel('Clustering coefficient - binary')
    ylabel('Count')
    title(sprintf('Clustering coefficient binary mean: %f; SD: %f', cc_bin, cc_bin_std))

    figure
    hist(local_efficiency_clean,10)
    xlabel('Local efficiency')
    ylabel('Count')
    title(sprintf('Local efficiency mean: %f; SD: %f', local_efficiency_mean, local_efficiency_std))
end
%% 7. Prepare output

clustering_parameter.cc.wei.collection=cc_values_clean;
clustering_parameter.cc.wei.mean=cc_mean_weighted;
clustering_parameter.cc.wei.std=cc_std_weighted;
clustering_parameter.cc.wei.distribution.x=cc_wei_x;
clustering_parameter.cc.wei.distribution.n=cc_wei_n;
clustering_parameter.cc.wei.distribution.h=cc_wei_h;

clustering_parameter.cc.bin.collection=cc_values_bin;
clustering_parameter.cc.bin.mean=cc_bin;
clustering_parameter.cc.bin.std=cc_bin_std;
clustering_parameter.cc.bin.distribution.x=cc_bin_x;
clustering_parameter.cc.bin.distribution.n=cc_bin_n;
clustering_parameter.cc.bin.distribution.h=cc_bin_h;

clustering_parameter.local_efficiency.collection=local_efficiency_clean;
clustering_parameter.local_efficiency.mean=local_efficiency_mean;
clustering_parameter.local_efficiency.std=local_efficiency_std;
clustering_parameter.local_efficiency.distribution.x=le_x;
clustering_parameter.local_efficiency.distribution.n=le_n;
clustering_parameter.local_efficiency.distribution.h=le_h;

clustering_parameter.transitivity=transitivity;

clustering_parameter.modularity_degree=modularity_degree;




end