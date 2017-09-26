%% Statistical analysis of network parameter as produced by runing analysis scripts A through D



data_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % data path to where to find the results of previous parts of the analysis
data_prefix='connect_sulfa_'; % filename stem of the results files


n_cond=2; % number of different experimental conditions
condition_labels={'CTRL', 'SAS'}; % names of the different experimental conditions
recordingnumbers.condition1=[1,3,6, 8,9,   13,14,16,17,20]; % which recordings are of experimental condition 1
recordingnumbers.condition2=[2,4,5,7, 10,11,12,  15,18,19,21,22]; % which recordings are of experimental condition 2
% add here more lines of recording conditions, respectively





%% Basic calculations

cell_parameter={'in_degree', 'out_degree', 'opi_degree', 'spike_rate'};
connection_parameter={'pp', 'pd'};
network_parameter={'netw_strength', 'connectivity_degree', 'cpl_bin', 'cpl_wei', 'global_efficiency', 'cc_bin', 'cc_wei', 'transitivity', 'modularity_degree', 'centrality', 'assortativity', 'NWF', 'kappa'};




% Check completeness of input
if n_cond ~=size(condition_labels,2)
    error('Number of conditions (n_cond) and number of condition labels (condition_labels) don''t match')
end
recordings_to_load=[];
for c=1:n_cond
    eval(sprintf('n_recordings_cond%i=size(recordingnumbers.condition%i,2);',c, c));  % if this line throws an error, check if recording numbers are specified correctly for all experimental conditions
    eval(sprintf('recordings_to_load=[recordings_to_load recordingnumbers.condition%i];', c));
end

% Initialize variables
for param_n=1:size(cell_parameter,2)
    param=cell_parameter{1,param_n};
    eval(sprintf('celldata.%s=[];', param));
end
celldata.recording=[];
celldata.condition=[];
for param_n=1:size(connection_parameter,2)
    param=connection_parameter{1,param_n};
    eval(sprintf('conndata.%s=[];', param));
end
conndata.recording=[];
conndata.condition=[];
for param_n=1:size(network_parameter,2)
    param=network_parameter{1,param_n};
    eval(sprintf('netwdata.%s=[];', param));
end
netwdata.recording=[];
netwdata.condition=[];

% Arrange figures
spacing=0.02;
table_height=0.3;
figure_height=0.8;
table_width=0.41+0.0446*n_cond;
figure_width=1-spacing-spacing-spacing-table_width;
table_position=[spacing, 0.4, table_width, table_height];
figure_position=[(spacing+spacing+table_width), (1-spacing-figure_height), figure_width, figure_height];
space_sum=19+2.1*n_cond;

%% Load data

n_recordings = max(structfun(@(x) max(x(:)), recordingnumbers));

for i=1:n_recordings
    if any(ismember(i, recordings_to_load))
        recordings_networks(i,1)=load(sprintf('%s%s%i_network_parameters.mat', data_path, data_prefix, i));
        recordings_activity(i,1)=load(sprintf('%s%s%i_network_activity.mat', data_path, data_prefix, i));        
        cond=find(structfun(@(x)any(ismember(x,i)),recordingnumbers));
        celldata.in_degree=[celldata.in_degree; recordings_networks(i,1).degree_parameter.in_degree.norm];
        celldata.out_degree=[celldata.out_degree; recordings_networks(i,1).degree_parameter.out_degree.norm];
        celldata.opi_degree=[celldata.opi_degree; recordings_networks(i,1).degree_parameter.opi_degree.collection;];
        celldata.spike_rate=[celldata.spike_rate; recordings_activity(i,1).spike_frequency];        
        celldata.recording=[celldata.recording; repmat(i,(size(celldata.in_degree,1)-size(celldata.recording,1)),1)];
        celldata.condition=[celldata.condition; repmat(condition_labels(1,cond),(size(celldata.in_degree,1)-size(celldata.condition,1)),1)];
        
        conndata.pp=[conndata.pp; recordings_networks(i,1).pp_parameter.pp.collection];
        conndata.pd=[conndata.pd; recordings_networks(i,1).pp_parameter.pd.collection];
        conndata.recording=[conndata.recording; repmat(i,(size(conndata.pp,1)-size(conndata.recording,1)),1)];
        conndata.condition=[conndata.condition; repmat(condition_labels(1,cond),(size(conndata.pp,1)-size(conndata.condition,1)),1)];
        
        netwdata.netw_strength=[netwdata.netw_strength; recordings_networks(i,1).pp_parameter.pp.mean];
        netwdata.connectivity_degree=[netwdata.connectivity_degree;  recordings_networks(i,1).pp_parameter.connectivity_degree];
        netwdata.cpl_bin=[netwdata.cpl_bin; recordings_networks(i,1).path_length_parameter.bin.cpl];
        netwdata.cpl_wei=[netwdata.cpl_wei; recordings_networks(i,1).path_length_parameter.wei.cpl];
        netwdata.global_efficiency=[netwdata.global_efficiency; recordings_networks(i,1).path_length_parameter.global_efficiency];
        netwdata.cc_bin=[netwdata.cc_bin; recordings_networks(i,1).clustering_parameter.cc.bin.mean];
        netwdata.cc_wei=[netwdata.cc_wei; recordings_networks(i,1).clustering_parameter.cc.wei.mean];
        netwdata.transitivity=[netwdata.transitivity; recordings_networks(i,1).clustering_parameter.transitivity];
        netwdata.modularity_degree=[netwdata.modularity_degree; recordings_networks(i,1).clustering_parameter.modularity_degree];
        netwdata.centrality=[netwdata.centrality; recordings_networks(i,1).structure_parameter.betweenness_centrality.mean];
        netwdata.assortativity=[netwdata.assortativity; recordings_networks(i,1).structure_parameter.assortativity];
        netwdata.NWF=[netwdata.NWF; recordings_activity(i,1).culture_wide_firing_rate];
        netwdata.kappa=[netwdata.kappa; recordings_activity(i,1).mean_kappa];        
        netwdata.recording=[netwdata.recording; repmat(i,(size(netwdata.netw_strength,1)-size(netwdata.recording,1)),1)];
        netwdata.condition=[netwdata.condition; repmat(condition_labels(1,cond),(size(netwdata.netw_strength,1)-size(netwdata.condition,1)),1)];     
    end
end

%% Descriptive statistics: cell parameter

for param_type=1:3
    switch param_type
        case (1)
            parameter_to_use=cell_parameter;
            data_to_use='celldata';
        case (2)
            parameter_to_use=connection_parameter;
            data_to_use='conndata';
        case (3)
            parameter_to_use=network_parameter;
            data_to_use='netwdata';
    end
    for param_n=1:size(parameter_to_use,2)
        param=parameter_to_use{1,param_n};
        for c=1:n_cond
            eval(sprintf('data=%s.%s(cellfun(@strcmp, %s.condition, repmat(condition_labels(1,c), size(%s.condition,1),1)));', data_to_use, param, data_to_use, data_to_use));
            eval(sprintf('statistics.%s.desc(3,c)=median(data);', param));
            eval(sprintf('statistics.%s.desc(6,c)=mean(data);', param));
            eval(sprintf('statistics.%s.desc(7,c)=std(data);', param));
            eval(sprintf('statistics.%s.desc([1,2,4,5],c)=prctile(data,[10;25;75;90]);', param));
            minimum=min(data);
            maximum=max(data);
            eval(sprintf('statistics.%s.desc(8,c)=minimum;',param));
            eval(sprintf('statistics.%s.desc(9,c)=maximum;', param));
            eval(sprintf('statistics.%s.desc(10,c)=range(data);', param));
            eval(sprintf('statistics.%s.desc(11,c)=skewness(data);', param));
            n=size(data,1);
            eval(sprintf('statistics.%s.desc(12,c)=n;', param));
            if n>10
                eval(sprintf('bin_size=2*(statistics.%s.desc(4,c)-statistics.%s.desc(2,c))/(n^(1/3));', param, param));
                edges=[minimum:bin_size:maximum];
                if size(edges,2)<16;
                    edges=linspace(min(data), max(data), 16);
                end
                [y,edges]=histcounts(data, edges);
                x=(edges(2:end)+edges(1:end-1))./2;
                try
                [fitobj, gof]=fit(x',y','gauss1', 'Lower', [-Inf,minimum, -Inf], 'Upper', [Inf, maximum, Inf], 'StartPoint', [rand(1,1), x(find(y(1,:)==max(y), 1)), ((x(1,end)-x(1,1))/2)]);
                eval(sprintf('statistics.%s.desc(13,c)=gof.rsquare;', param));
                catch
                    eval(sprintf('statistics.%s.desc(13,c)=0;', param));
                    warning(sprintf('Normal fit for %s %s failed, returning 0 as R2', condition_labels{1,c}, param))
                end
            else
                eval(sprintf('statistics.%s.desc(13,c)=0;', param));
                warning('Goodness of normal fit can only be assed for more than 9 samples')
            end
        end
    end
end
table_row_names={'10% percentile', '25% percentile', 'median', '75% percentile', '90% percentile', 'mean', 'std', 'min', 'max', 'range', 'skewness', '# of samples', 'R2 normal fit'};

%% Statistic inference: cell parameter
for param_n=1:size(cell_parameter,2)
    param=cell_parameter{1,param_n};
    comparison_n=1;
    comparisons_to_make=triu(ones(n_cond,n_cond),1);
    comparison_labels={};
    for c1=1:n_cond
        for c2=1:n_cond
            if comparisons_to_make(c1,c2)==1
                eval(sprintf('data1=celldata.%s(cellfun(@strcmp, celldata.condition, repmat(condition_labels(1,c1), size(celldata.condition,1),1)));', param));
                eval(sprintf('data2=celldata.%s(cellfun(@strcmp, celldata.condition, repmat(condition_labels(1,c2), size(celldata.condition,1),1)));', param));
                mes_obj=mes(data1, data2, 'md');
                eval(sprintf('statistics.%s.comp{comparison_n,1}=mes_obj.md;', param));
                mes_obj=mes(data1, data2, 'U1');
                eval(sprintf('statistics.%s.comp{comparison_n,2}=mes_obj.U1;', param));
                eval(sprintf('statistics.%s.comp{comparison_n,3}=ranksum(data1, data2);', param));
                eval(sprintf('statistics.%s.comp{comparison_n,4}=''Wilcoxon'';', param));
                comparison_labels(1,comparison_n)={sprintf('%s vs. %s', condition_labels{1,c1}, condition_labels{1,c2})};
                comparison_n=comparison_n+1;
            end
        end
    end
    
end
table_column_names_cell={'mean difference', 'Cohens U', 'p', 'test'};

%% Statistic inference: connection parameter
for param_n=1:size(connection_parameter,2)
    param=connection_parameter{1,param_n};
    comparison_n=1;
    comparisons_to_make=triu(ones(n_cond,n_cond),1);
    comparison_labels={};
    for c1=1:n_cond
        for c2=1:n_cond
            if comparisons_to_make(c1,c2)==1
                eval(sprintf('data1=conndata.%s(cellfun(@strcmp, conndata.condition, repmat(condition_labels(1,c1), size(conndata.condition,1),1)));', param));
                eval(sprintf('data2=conndata.%s(cellfun(@strcmp, conndata.condition, repmat(condition_labels(1,c2), size(conndata.condition,1),1)));', param));
                mes_obj=mes(data1, data2, 'md');
                eval(sprintf('statistics.%s.comp{comparison_n,1}=mes_obj.md;', param));
                mes_obj=mes(data1, data2, 'U1');
                eval(sprintf('statistics.%s.comp{comparison_n,2}=mes_obj.U1;', param));
                eval(sprintf('statistics.%s.comp{comparison_n,3}=ranksum(data1, data2);', param));
                eval(sprintf('statistics.%s.comp{comparison_n,4}=''Wilcoxon'';', param));
                comparison_labels(1,comparison_n)={sprintf('%s vs. %s', condition_labels{1,c1}, condition_labels{1,c2})};
                comparison_n=comparison_n+1;
            end
        end
    end
    
end
table_column_names_conn={'mean difference', 'Cohens U', 'p', 'test'};

%% Statistic inference: network parameter
for param_n=1:size(network_parameter,2)
    param=network_parameter{1,param_n};
    comparison_n=1;
    comparisons_to_make=triu(ones(n_cond,n_cond),1);
    comparison_labels={};
    for c1=1:n_cond
        for c2=1:n_cond
            if comparisons_to_make(c1,c2)==1
                eval(sprintf('data1=netwdata.%s(cellfun(@strcmp, netwdata.condition, repmat(condition_labels(1,c1), size(netwdata.condition,1),1)));', param));
                eval(sprintf('data2=netwdata.%s(cellfun(@strcmp, netwdata.condition, repmat(condition_labels(1,c2), size(netwdata.condition,1),1)));', param));
                mes_obj=mes(data1, data2, 'md');
                eval(sprintf('statistics.%s.comp{comparison_n,1}=mes_obj.md;', param));
                mes_obj=mes(data1, data2, 'U1');
                eval(sprintf('statistics.%s.comp{comparison_n,2}=mes_obj.U1;', param));
                mes_obj=mes(data1,data2, 'glassdelta');
                eval(sprintf('statistics.%s.comp{comparison_n,3}=mes_obj.glassdelta;', param));
                [test_obj(1,1), test_obj(1,2)]=ttest2(data1, data2);
                eval(sprintf('statistics.%s.comp{comparison_n,4}=test_obj(1,2);', param));
                eval(sprintf('statistics.%s.comp{comparison_n,5}=''t-test'';', param));
                comparison_labels(1,comparison_n)={sprintf('%s vs. %s', condition_labels{1,c1}, condition_labels{1,c2})};
                comparison_n=comparison_n+1;
            end
        end
    end
    
end
table_column_names_netw={'mean difference', 'Cohens U', 'Glass delta', 'p', 'test'};

%% Plot cell parameter
for param_n=1:size(cell_parameter,2)
    clear g
    param=cell_parameter{1,param_n};
    
    eval(sprintf('figure_tables = figure(''Units'', ''normalized'', ''Position'', table_position, ''name'', ''statistics %s'');', param));
    eval(sprintf('table_desc = uitable(figure_tables,''Data'',statistics.%s.desc,''ColumnName'',condition_labels,''RowName'', table_row_names ,''Units'', ''normalized'', ''Position'',[0, 0, (4.3+2.1*n_cond)/space_sum, 1]);', param));
    eval(sprintf('table_comp = uitable(figure_tables,''Data'',statistics.%s.comp,''ColumnName'',table_column_names_cell,''RowName'', comparison_labels ,''Units'', ''normalized'', ''Position'',[(4.3+2.1*n_cond)/space_sum, 0, (14.7)/space_sum, 1]);', param));
    
    eval(sprintf('g(2,1)=gramm(''x'', celldata.%s);', param));
    g(2,1).facet_grid(celldata.condition,[]);
    g(2,1).set_continuous_color('active',false);
    g(2,1).stat_bin();
    g(2,1).set_title('All values');
    g(2,1).set_order_options('row',condition_labels);      
    eval(sprintf('g(2,1).set_names(''x'', ''%s'', ''color'', ''rec'', ''row'', '''', ''y'', '''');', param));
    
    eval(sprintf('g(1,1)=gramm(''x'', celldata.%s, ''color'', celldata.recording);', param));
    g(1,1).facet_grid(celldata.condition,[]);
    g(1,1).stat_density('kernel', 'triangle');
    g(1,1).set_title('Single recordings');
    g(1,1).set_continuous_color('active',false);
    g(1,1).set_order_options('row',condition_labels);      
    eval(sprintf('g(1,1).set_names(''x'', ''%s'', ''color'', ''rec'', ''row'', '''', ''y'', '''');', param));
    
    eval(sprintf('g(2,2)=gramm(''x'', celldata.condition, ''y'', celldata.%s);', param));
    g(2,2).stat_boxplot('width',0.15);
    g(2,2).set_title('All values');
    g(2,2).set_order_options('x',condition_labels);      
    eval(sprintf('g(2,2).set_names(''x'', ''condition'', ''y'', ''%s'');', param));
    
    eval(sprintf('g(1,2)=gramm(''x'', celldata.condition, ''y'', celldata.%s, ''color'', celldata.recording);', param));
    g(1,2).geom_jitter();
    g(1,2).set_title('Single values');
    g(1,2).set_continuous_color('active',false);
    g(1,2).set_order_options('x',condition_labels);      
    eval(sprintf('g(1,2).set_names(''x'', ''condition'', ''y'', ''%s'', ''color'', ''rec'');', param));
    
    eval(sprintf('figure(''units'',''normalized'',''outerposition'', figure_position, ''name'', ''Visualization %s'');', param));
    g.draw();
    
end

%% Plot connection parameter
for param_n=1:size(connection_parameter,2)
    clear g
    param=connection_parameter{1,param_n};
    eval(sprintf('figure_tables = figure(''Units'', ''normalized'', ''Position'', table_position, ''name'', ''statistics %s'');', param));
    eval(sprintf('table_desc = uitable(figure_tables,''Data'',statistics.%s.desc,''ColumnName'',condition_labels,''RowName'', table_row_names ,''Units'', ''normalized'', ''Position'',[0, 0, (4.3+2.1*n_cond)/space_sum, 1]);', param));
    eval(sprintf('table_comp = uitable(figure_tables,''Data'',statistics.%s.comp,''ColumnName'',table_column_names_conn,''RowName'', comparison_labels ,''Units'', ''normalized'', ''Position'',[(4.3+2.1*n_cond)/space_sum, 0, (14.7)/space_sum, 1]);', param));
    
    eval(sprintf('g(2,1)=gramm(''x'', conndata.%s);', param));
    g(2,1).facet_grid(conndata.condition,[]);
    g(2,1).set_continuous_color('active',false);
    g(2,1).stat_bin();
    g(2,1).set_title('All values');
    g(2,1).set_order_options('row',condition_labels);      
    eval(sprintf('g(2,1).set_names(''x'', ''%s'', ''color'', ''rec'', ''row'', '''', ''y'', '''');', param));
    
    eval(sprintf('g(1,1)=gramm(''x'', conndata.%s, ''color'', conndata.recording);', param));
    g(1,1).facet_grid(conndata.condition,[]);
    g(1,1).stat_density('kernel', 'triangle');
    g(1,1).set_title('Single recordings');
    g(1,1).set_continuous_color('active',false);
    g(1,1).set_order_options('row',condition_labels);      
    eval(sprintf('g(1,1).set_names(''x'', ''%s'', ''color'', ''rec'', ''row'', '''', ''y'', '''');', param));
    
    eval(sprintf('g(2,2)=gramm(''x'', conndata.condition, ''y'', conndata.%s);', param));
    g(2,2).stat_boxplot('width',0.15);
    g(2,2).set_title('All values');
    g(2,2).set_order_options('x',condition_labels);      
    eval(sprintf('g(2,2).set_names(''x'', ''condition'', ''y'', ''%s'');', param));
    
    eval(sprintf('g(1,2)=gramm(''x'', conndata.condition, ''y'', conndata.%s, ''color'', conndata.recording);', param));
    g(1,2).geom_jitter();
    g(1,2).set_title('Single values');
    g(1,2).set_continuous_color('active',false);
    g(1,2).set_order_options('x',condition_labels);      
    eval(sprintf('g(1,2).set_names(''x'', ''condition'', ''y'', ''%s'', ''color'', ''rec'');', param));
    
    eval(sprintf('figure(''units'',''normalized'',''outerposition'', figure_position, ''name'', ''Visualization %s'');', param));
    g.draw();
    
end

%% Plot network parameter

for param_n=1:size(network_parameter,2)
    clear g
    param=network_parameter{1,param_n};
    eval(sprintf('figure_tables = figure(''Units'', ''normalized'', ''Position'', table_position, ''name'', ''statistics %s'');', param));
    eval(sprintf('table_desc = uitable(figure_tables,''Data'',statistics.%s.desc,''ColumnName'',condition_labels,''RowName'', table_row_names ,''Units'', ''normalized'', ''Position'',[0, 0, (4.3+2.1*n_cond)/space_sum, 1]);', param));
    eval(sprintf('table_comp = uitable(figure_tables,''Data'',statistics.%s.comp,''ColumnName'',table_column_names_netw,''RowName'', comparison_labels ,''Units'', ''normalized'', ''Position'',[(4.3+2.1*n_cond)/space_sum, 0, (14.7)/space_sum, 1]);', param));
    
    eval(sprintf('g(1,1)=gramm(''x'', netwdata.condition, ''y'', netwdata.%s, ''color'', netwdata.recording);', param));
    g(1,1).geom_jitter();
    g(1,1).set_title('Single values');
    g(1,1).set_continuous_color('active',false);
    g(1,1).set_order_options('x',condition_labels);      
    eval(sprintf('g(1,1).set_names(''x'', ''condition'', ''y'', ''%s'', ''color'', ''rec'');', param));
    
    eval(sprintf('g(1,2)=gramm(''x'', netwdata.condition, ''y'', netwdata.%s);', param));
    g(1,2).stat_violin('normalization','area','dodge',0,'fill','edge');
    g(1,2).stat_boxplot('width',0.15);
    g(1,2).set_title('All values');
    g(1,2).set_order_options('x',condition_labels);      
    eval(sprintf('g(1,2).set_names(''x'', ''condition'', ''y'', ''%s'');', param));
    
    eval(sprintf('figure(''units'',''normalized'',''outerposition'', figure_position, ''name'', ''Visualization %s'');', param));
    g.draw();
    
end
