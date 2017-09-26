close all
clearvars -except TrainedClassifier

data_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % data path to find the results of part C
data_filename='connect_sulfa'; % filename stem of results of part C
objective=4; % which objective 10x, 4x, 2x, ... (Note that conversion from image pixels to µm might differ from the used factors with your individual microscopy setup)
show_figures=0; % Show figures?
max_percent_connected=100; % Change this if you want to analyze only the stronges x% links in the network (This is not generally recommended)
savedata=1; % Save results?
save_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % Data path to where to save the results
do_all=1;  % quantify all network topology parameters?
do='degree'; % quatify a specific set of network topology parameters (only relevant if do_all=0)
% degree=In-degree (number of incomming connections per cell); out-degree (number of outgoing connections per cell); out-per-in-degree (ratio of outgoing per incomming connections per cell)
% pp_pd= propagation propability and physical distance
% pl=path lengths
% cluster= clustering coefficient, transitivity and local efficiency
% structure= analysis on betweenness centrality, assortativity coefficient, small-worldness

messungen=[13:1:22]; % batch of recordings that is to be analyzed



%% 0. complete inputs
set(0, 'DefaulttextInterpreter', 'none')
dispstat('', 'init')
number_of_scripts=5;
todo=zeros(number_of_scripts,1);
if do_all==1
    clearvars do
    todo(:,1)=ones(number_of_scripts,1);
else
    if strcmp(do, 'degree')==1
        todo(1,1)=1;
    elseif strcmp(do, 'pp_pd')==1
        todo(2,1)=1;
    elseif strcmp(do, 'pl')==1
        todo(3,1)=1;
    elseif strcmp(do, 'cluster')==1
        todo(4,1)=1;
    elseif strcmp(do, 'structure')==1
        todo(3,1)=1;
        todo(4,1)=1;
        todo(5,1)=1;
    end
end

micrometer_per_pixel=22.8/objective;
dispstat('Loading classification model', 'timestamp')
TrainedClassifier=load('Trained_classifier_final8.mat'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% If relevant imaging parameters differ in your setup you must retrain the classification model on simulations adopted to your individual imaging parameters


%%
for counter=1:size(messungen,2)
    %     close all
    dispstat(sprintf('Start network analysis for measurement %i', messungen(counter)), 'timestamp', 'keepthis')
    name=sprintf('%s_%i',data_filename, messungen(counter));
    
    %% 1. Load data
    dispstat('Loading data', 'timestamp')
    corr_name_norm=sprintf('%s_reconstructed_network_norm.mat', name);
    correlation_values_norm=load([data_path corr_name_norm]);
    
    corr_name=sprintf('%s_reconstructed_network.mat', name);
    correlation_values=load([data_path corr_name]);
    
    reg_name=sprintf('%s_region_properties.mat', name);
    load([data_path reg_name])
    n_n=size(region_properties,1);

    
    

    
    
    
    
    
    
    
    %% 2. Arrange predictors
    dispstat('Preparing predictors', 'timestamp')
    [ predictors ] = arrange_predictors(  correlation_values_norm.network_xcorr, correlation_values_norm.network_MI, correlation_values_norm.network_JE, correlation_values_norm.network_TE, correlation_values_norm.network_GTE, correlation_values_norm.network_SC, correlation_values_norm.network_PP, correlation_values_norm.network_PP_samebin );
    [ predictors_nonorm ] = arrange_predictors(  correlation_values.network_xcorr, correlation_values.network_MI, correlation_values.network_JE, correlation_values.network_TE, correlation_values.network_GTE, correlation_values.network_SC, correlation_values.network_PP, correlation_values.network_PP_samebin );
    
    
    %% 3. Apply model
    clear Yfit
    dispstat('Applying classifation model', 'timestamp')
    [Yfit(:,1), Yfit(:,2:3)] = predict(TrainedClassifier.rusTree,predictors);  % Yfit = predicted classes which we believe to be "true"
    
    counter2=1;
    for s=1:n_n
        for t=1:n_n
            if s~=t
                Yfit(counter2,4)=s;
                Yfit(counter2,5)=t;
                counter2=counter2+1;
            end
        end
    end
    Yfit(:,6)=predictors_nonorm(:,8);   
    
    Yfit=Yfit(Yfit(:,1)==1,:);  %only use pairs classified as connected
    
    
    %Rebulid connection matrix
    maxscore=max(Yfit(:,3));
    minscore=min(Yfit(:,3));
    if maxscore~=0
        scores_normalized=(Yfit(:,3)-minscore)./maxscore;
    else
        scores_normalized=Yfit(:,3)-minscore;
    end
    correlation_matrix=zeros(n_n,n_n);
    correlation_matrix_bin=zeros(n_n,n_n);
    
    for conn=1:size(Yfit,1)
        if scores_normalized(conn,1)>=(1-max_percent_connected/100)
            correlation_matrix(Yfit(conn,4),Yfit(conn,5))=Yfit(conn,6);
            correlation_matrix_bin(Yfit(conn,4),Yfit(conn,5))=1;
        end
        conn=conn+1;
    end
    
    

        save(sprintf('%s%s_network.mat', save_path, name), 'correlation_matrix')
        dispstat(sprintf('Network saved as %s_network.mat',name), 'timestamp', 'keepthis')

    
    
    dispstat('Classification complete', 'timestamp', 'keepthis')
    %% 4. Network analysis
    dispstat('Analyzing the network ...', 'timestamp', 'keepthis')
    %% 4.1. Basic calculations
    
    connection_length_matrix=1./correlation_matrix;
    [distance_matrix, distance_matrix_bin]=distance_wei(connection_length_matrix);
    % correlation_matrix = weighted network
    % correlation_matrix_bin = binary network
    % distance_matrix = distance_matrix of weighted network
    % distance_matrix_bin = binary distance matrix of binary network
    
    
    number_of_cells=size(correlation_values_norm.network_xcorr,1);
    
    for cell=1:number_of_cells
        coors(cell,1:2)=region_properties(cell,1).Centroid(1,1:2);
    end
    
    % Remove unconnected cells
    for i=1:number_of_cells
        region_propertiees(i,1).id=i;
    end
    cells_in_network=[1:1:number_of_cells]';    
    for cell=1:number_of_cells
        conn_sum(cell,1)=sum(correlation_matrix_bin(cell,:))+sum(correlation_matrix_bin(:,cell));
        connected_cells(cell,1)=cell;
    end
    for cell=number_of_cells:-1:1
        if sum(correlation_matrix_bin(:,cell))+sum(correlation_matrix_bin(cell,:))==0
            correlation_matrix_bin(cell,:)=[];
            correlation_matrix_bin(:,cell)=[];
            correlation_matrix(cell,:)=[];
            correlation_matrix(:,cell)=[];
            connected_cells(cell,:)=[];
                        coors(cell,:)=[];
            region_properties(cell,:)=[];
            cells_in_network(cell,:)=[];
            number_of_cells=number_of_cells-1;
        end
    end
    
    % Built colormap
    number_of_connections=sum(correlation_matrix_bin(:));
    cmap=ones(number_of_connections,3);
    cmap(:,1)=cmap(:,1).*0.2;
    cmap(:,3)=cmap(:,3).*0.4;
    cmap(:,2)=(Yfit(:,3)-min(Yfit(:,3)));
    cmap(:,2)=cmap(:,2)/max(cmap(:,2));
    
    % Connectivity degree
    real_connectivity_degree=number_of_connections/(number_of_cells*(number_of_cells-1));
    
    
        if show_figures==1
    figure('units','normalized','outerposition',[0 0 1 1])
    %         subplot(1,2,1)
    colormap('bone')
    imagesc(diff_image)
    hold on
    gplotcolor(correlation_matrix_bin,coors,cmap)
    axis ij
    hold on
    scatter(coors(:,1), coors(:,2),10, 'white')
            title(sprintf('%s network', name))
    pbaspect([1,1,1])
    
    
    figure
    % subplot(1,2,2)
    colormap('bone')
    imagesc(diff_image)
    pbaspect([1,1,1])
                end
    drawnow
    %% 4.2. Network analysis scripts
    
    
    
    % calculation of quantitative parameters
    if todo(1,1)==1
        [degree_parameter]=degree_analysis(number_of_cells, correlation_matrix_bin, region_properties, show_figures);
        dispstat('Degree analysis complete', 'timestamp', 'keepthis')
    end
    if todo(2,1)==1
        [pp_parameter]=pp_analysis(correlation_matrix, correlation_matrix_bin, number_of_connections, number_of_cells,real_connectivity_degree, micrometer_per_pixel, region_properties, show_figures);
        dispstat('Propagation analysis complete', 'timestamp', 'keepthis')
    end
    if todo(3,1)==1
        [path_length_parameter]=path_length_calculations_vul(distance_matrix_bin, distance_matrix, correlation_matrix, show_figures);
        dispstat('Path length analysis complete', 'timestamp', 'keepthis')
    end
    if todo(4,1)==1
        [clustering_parameter]=clustering_analysis(correlation_matrix, number_of_connections, correlation_matrix_bin, show_figures);
        dispstat('Clustering analysis complete', 'timestamp', 'keepthis')
    end
    if todo(5,1)==1
        [structure_parameter]=structure_analysis(distance_matrix, number_of_connections, correlation_matrix, show_figures);
        dispstat('Structure analysis complete', 'timestamp', 'keepthis')
    end
    
    
    %% 5. Save results
    % savemultfigs(name, algorithm)
    
    if savedata==1
        save(sprintf('%s%s_network_parameters.mat', save_path, name), 'degree_parameter', 'pp_parameter', 'path_length_parameter', 'clustering_parameter', 'structure_parameter')
        dispstat(sprintf('Data saved as %s_network_parameters.mat',name), 'timestamp', 'keepthis')
        
        save(sprintf('%s%s_region_properties_network.mat', save_path, name), 'region_properties', 'cells_in_network')
        dispstat(sprintf('Region properties saved as %s_region_properties_network.mat', name), 'timestamp', 'keepthis')        
        
    end
    
    
    
    clearvars -except connectivity_degree counter data_path do do_all messungen name objective save_path savedata show_figures simulationen n_n TrainedClassifier max_percent_connected todo micrometer_per_pixel
    
    
    dispstat(sprintf('Finished network analysis for measurement %i\n\n',messungen(counter)), 'timestamp', 'keepthis')
% pause
end

