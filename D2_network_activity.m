% close all

%% Input
rec_nr_start=13; % First recording to analyze when running batches
rec_nr_stop=22; % Last recording to analyze when running batches
data_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % Data path to where to find the results of parts C and D
data_filename='connect_sulfa'; % filename stem for the results of parts C and D
show_figures=0; % Show figures?
clear_data=1; %clear simulation data after each run (to save memory); will always be on if more than 20 simulations are batch processed
sampling_frequency=27.33; % sampling frequency of the recording in Hz
seconds_per_bin=1; %how many seconds should be summarized into one bin for the images



%% Initialize
data=[];
if rec_nr_stop-rec_nr_start>20
    clear_data=1;
    show_figures=0;
end



number_of_frames_per_bin=ceil(sampling_frequency*seconds_per_bin);

dispstat('','init')
%% Calculations
for rec_nr=rec_nr_start:rec_nr_stop
    dispstat(sprintf('Loading data for recording %i',rec_nr), 'keepthis', 'timestamp')
    filename_spikes=sprintf('%s%s_%i_estimated_spikes.mat', data_filename,data_path, rec_nr);
    filename_regions=sprintf('%s%s_%i_region_properties_network.mat', data_filename,data_path, rec_nr);    
    
    spikes_temp=load(filename_spikes);
    data(rec_nr,1).estimated_spikes=spikes_temp;
    regions_temp=load(filename_regions);
    data(rec_nr,1).cells_in_network=regions_temp.cells_in_network;   
    
    % Remove spike data from cells that are not part of the network

    data(rec_nr,1).estimated_spikes.raster_plot=data(rec_nr,1).estimated_spikes.raster_plot(data(rec_nr,1).cells_in_network,:);
    data(rec_nr,1).estimated_spikes.spike_trains=data(rec_nr,1).estimated_spikes.spike_trains(data(rec_nr,1).cells_in_network,:);
    
    
    
    n_cells=size(data(rec_nr,1).estimated_spikes.raster_plot,1);

    
    dispstat('Data passed to network activity analysis', 'keepthis', 'timestamp')
    %% Network activit analysis
    
    % bin data
    data_to_bin=data(rec_nr,1).estimated_spikes.raster_plot;
    data_binned=zeros(size(data_to_bin,1),1);
    
    for frame=1:number_of_frames_per_bin:size(data_to_bin,2)-number_of_frames_per_bin
        data_binned(:,end+1)=sum(data_to_bin(:,frame:frame+number_of_frames_per_bin-1),2);
    end
        n_bins=size(data_binned,2);
    
    
    %% Generate raster plot
    
    if show_figures==1
        figure('units','normalized','outerposition',[0 0 1 1])
        s1=subplot(5,1,(1:4));
        imagesc(data_binned);
        set(gca, 'YGrid', 'on')
        ybounds = ylim;
        set(gca,'YTick',ybounds(1):ybounds(2));
        labels=get(gca, 'YTickLabel');
        for tick=2:size(labels,1)-1
            labels{tick,1}='';
        end
        set(gca, 'YTickLabel', labels, 'XTickLabel', [])
        title('data binned')
        ylabel('cells')
        
        s2=subplot(5,1,5);
        plot(sum(data_binned,1), 'Color', 'black');
        linkaxes([s1,s2], 'x')
        cmap=colormap('bone');
        colormap(flipud(cmap));
        ylabel('spikes/bin')
        xlabel('bins')
        drawnow
    end
    
    %% Calculate spike frequency / firing rate
    
    recording_duration=size(data_to_bin,2)/sampling_frequency;
    spike_frequency=sum(data(rec_nr,1).estimated_spikes.raster_plot,2)/recording_duration; %spike frequency in Hz
    spike_frequency_frames=sum(data(rec_nr,1).estimated_spikes.raster_plot,2)/size(data(rec_nr,1).estimated_spikes.raster_plot,2); %spike frequency in 1/frame
    culture_wide_firing_rate=mean(spike_frequency);
       
    
    
    %% Calculate burst frequency https://bmcneurosci.biomedcentral.com/articles/10.1186/1471-2202-10-93
    for cell=1:size(data(rec_nr,1).estimated_spikes.raster_plot,1)
        data(rec_nr,1).estimated_spikes.burst_times{cell,1}(1,1)=0;
        if size(data(rec_nr,1).estimated_spikes.spike_trains{cell,1},1)>3
            for burst=1:size(data(rec_nr,1).estimated_spikes.spike_trains{cell,1},1)-2
                delta_frames=data(rec_nr,1).estimated_spikes.spike_trains{cell,1}(burst+2,1)-data(rec_nr,1).estimated_spikes.spike_trains{cell,1}(burst,1);
                p=poisspdf(3,spike_frequency_frames(cell,1)*delta_frames);
                if p<0.05
                    burst_times{cell,1}=[data(rec_nr,1).estimated_spikes.burst_times{cell,1}; data(rec_nr,1).estimated_spikes.spike_trains{cell,1}(burst,1)];
                end
            end
        end
    end
    
    
    
    
    %% Calculate measure of synchrony (Cohen's kappa) 
    
    
    comparisons_to_make=triu(ones(n_cells,n_cells),1);   
    kappa=zeros(size(comparisons_to_make));
    for cell1=1:n_cells
        s1=sum(data_binned(cell1,:),2);
        for cell2=1:n_cells
            if comparisons_to_make(cell1,cell2)==1
                        s2=sum(data_binned(cell2,:),2);
                        runter_summe=sum(data_binned([cell1, cell2],:),1);
                        c12=sum( runter_summe==2);
                        p_obs=(2*c12+n_bins-s1-s2)/(n_bins);
                        p_exp=(((s1*s2)/(n_bins))+(((n_bins-s1)*(n_bins-s2))/(n_bins)))/(n_bins);
                        kappa(cell1,cell2)=(p_obs-p_exp)/(1-p_exp);
            end
        end
    end
                
    kappa=kappa(logical(comparisons_to_make));    
    
    mean_kappa=mean(kappa);
            
    
    

    
    
    
    
    
    
    
        
    filename_save=sprintf('%sdata_filename_%i_network_activity.mat', data_filename,data_path, rec_nr);
        
    save(filename_save, 'spike_frequency', 'culture_wide_firing_rate', 'burst_times', 'kappa', 'mean_kappa');
    
    
    
    
    
    if clear_data==1
        data(rec_nr,:)=[];
    end
    
    
    dispstat(sprintf('Network activity analyzed for recording %i (Currently running recording %i through %i)\n', rec_nr, rec_nr_start, rec_nr_stop), 'keepthis', 'timestamp')
    
    
end





