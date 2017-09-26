close all


rec_nr_start=13; % First recording to be analyzed when running batches
rec_nr_stop=22; % Last recording to be analyzed when running batches
method='threshold2'; % Keep this input
data_path='V:\AG_Neurophotonik\Projekte\CONNECT-Sulfasalazine\Data\Networks\'; % Data path to find the results of part A
show_figures=0; % Show figures?
clear_data=1; %clear simulation data after each run (to save memory); will always be on if more than 20 simulations are batch processed
data_filename='connect_sulfa'; % filename_stem of results of part A


%% _________

data=[];
if rec_nr_stop-rec_nr_start>20
    clear_data=1;
    show_figures=0;
end




dispstat('','init')
%%
for rec_nr=rec_nr_start:rec_nr_stop
    dispstat(sprintf('Loading data for recording %i',rec_nr), 'keepthis', 'timestamp')
    filename_meanstack=sprintf('%s%s_%i_meanstack.mat', data_filename, data_path, rec_nr);

    meanstack_temp=load(filename_meanstack);
    data(rec_nr,1).meanstack=meanstack_temp.meanstack;
    
    
    %% Prepare and rearrange data
    
    
    [data(rec_nr,1).dF_F]=relative_fluorescence(data(rec_nr,1).meanstack);
    dF_F=data(rec_nr,1).dF_F;
    save(sprintf('%s%s_%i_relative_fluorescence.mat', data_filename, data_path, rec_nr), 'dF_F');
    
    
    if show_figures==1
        for i=1:size(dF_F,1)
            figure;plot(dF_F(i,:))
            title(sprintf('Cell %i', i))
        end
    end
    
    clearvars dF_F
    
    % %% Temporary insert to use the already saved relative fluorescence
    %     fl_temp=load(sprintf('V:\\AG_Neurophotonik\\Projekte\\disuse_hypersensitivity\\Daten_disuse_ML\\disuse_%i_relative_fluorescence.mat', rec_nr));
    %     data(rec_nr,1).dF_F=fl_temp.dF_F;
    
    %% Generation of fluorescence traces and spike estimation on them
    n_n=size(data(rec_nr,1).dF_F,1);
    n_frames=size(data(rec_nr,1).dF_F,2);
    
     
    if strcmp(method, 'threshold2')==1
        [data(rec_nr,1).estimated_spikes] = spike_estimation_th2( data(rec_nr,1).dF_F(:,200:29000), show_figures );
        spike_trains=data(rec_nr,1).estimated_spikes.spike_trains;
        raster_plot=data(rec_nr,1).estimated_spikes.raster_plot;
        save(sprintf('%s%s_%i_estimated_spikes.mat',data_filename, data_path, rec_nr), 'spike_trains', 'raster_plot', 'method');
    end
    
    if clear_data==1
        data(rec_nr,:)=[];
    end
    
    dispstat(sprintf('Spike estimation completed for recording %i (Currently running numbers %i through %i)\n', rec_nr, rec_nr_start, rec_nr_stop), 'keepthis', 'timestamp')
    
    clearvars spikes
    
    
    
end





