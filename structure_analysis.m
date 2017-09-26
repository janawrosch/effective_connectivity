function [structure_parameter]=structure_analysis(distance_matrix, cutoffzeile, correlation_matrix, show_figures)
%% 1. Betweenness centrality
betweenness_centrality=betweenness_wei(distance_matrix);
betweenness_centrality_clean=betweenness_centrality(betweenness_centrality>0);
if isempty(betweenness_centrality_clean)==1
    betweenness_centrality_clean=0;
end
betweenness_centrality_mean=mean(betweenness_centrality_clean);
betweenness_centrality_std=std(betweenness_centrality_clean);

[bc_n, bc_x]=hist(betweenness_centrality_clean, ceil(cutoffzeile/3));
bc_h=kstest(bc_n);

%% 2. Assortativity coefficient
assortativity_coefficient=assortativity_wei(correlation_matrix, 1);

% %% 3. Small-worldness
% 
% directed_ws_model
% 
% if show_figures==1
% figure
% [ax,h1,h2]=plotyy(b_values, ws_cpl_mean,b_values,ws_cc_mean, 'semilogx', 'semilogx');
% set(h1(1), 'Color', 'red');
% set(h2(1), 'Color', 'blue');
% set(ax(1), 'YColor', 'red');
% set(ax(2), 'YColor', 'blue');
% xlabel('Rewiring propability b')
% ylabel(ax(1), 'relative characteristic path legth L(p)/L(1)')
% ylabel(ax(2), 'relative clustering coefficient C(p)/C(1)')
% 
% title(sprintf('Small world chart\nL=%i; p_0=%.2f; number of realizations =%i\nMeasured values are indicated as horizontal lines', number_of_cells, p0,realization))
% legend({'Characteristic path length', 'Clustering coefficient'})
% cpl_line=ones(size(b_values)).*(cpl_bin/ws_cpl(22,1));
% cc_line=ones(size(b_values)).*(cc_bin/ws_cc(22,1));
% hold(ax(1))
% plot(ax(1), b_values, cpl_line, 'Color', 'red')
% hold(ax(2))
% plot(ax(2), b_values, cc_line, 'Color', 'blue')
% hold(ax(2))
% hold(ax(1))
% ylim(ax(2), [min(cc_line(1,1)-0.1, 1) max(ws_cc_mean(1,1)+0.1, cc_line(1,1)+0.1)])
% ylim(ax(1), [min(cpl_line(1,1)-0.1, 1) max(ws_cpl_mean(1,1)+0.1, cpl_line(1,1)+0.1)])
% end

%% 3. Show figures

if show_figures==1
    figure
    bar(bc_x, bc_n)
    xlabel('betweenness centrality')
    ylabel('count')
    title(sprintf('Betweenness Centrality mean: %f; SD: %f', betweenness_centrality_mean, betweenness_centrality_std))
end
%% 4. Prepare output

structure_parameter.betweenness_centrality.collection=betweenness_centrality_clean;
structure_parameter.betweenness_centrality.mean=betweenness_centrality_mean;
structure_parameter.betweenness_centrality.std=betweenness_centrality_std;
structure_parameter.betweenness_centrality.distribution.x=bc_x;
structure_parameter.betweenness_centrality.distribution.n=bc_n;
structure_parameter.betweenness_centrality.distribution.h=bc_h;

structure_parameter.assortativity=assortativity_coefficient;




end