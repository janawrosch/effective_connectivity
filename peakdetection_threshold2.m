function [  estimated_spikes ] = peakdetection_threshold2( traces_raw, noise, number_of_images, number_of_cells, sprungbreite, maxbreite, show_figures)

spike_times_binary=zeros(size(traces_raw));
spike_times_weighted=zeros(size(traces_raw));
spike_amplitude_log(number_of_cells,1).frame(1,1)=0;
spike_amplitude_log(number_of_cells,1).amplitude(1,1)=0;
for cell=1:number_of_cells
    spike_amplitude_log(cell,1).frame(1,1)=0;
    spike_amplitude_log(cell,1).amplitude(1,1)=0;
    
    trace=traces_raw(cell,:);
    spike_count=0;
    for frame=1+sprungbreite:number_of_images-maxbreite
        
        temp_sprunghoehe=max(trace(frame-sprungbreite:frame))-min(trace(frame-sprungbreite:frame));
        if temp_sprunghoehe>= noise(cell,1)
            ausschnitt=trace(frame-maxbreite:frame+maxbreite);
            [~,I]=max(ausschnitt);
            if I==maxbreite+1
                
                spike_count=spike_count+1;
                spike_times_binary(cell,frame)=1;
                spike_times_weighted(cell,frame)=temp_sprunghoehe;
                spike_amplitude_log(cell,1).frame(spike_count,1)=frame;
                spike_amplitude_log(cell,1).amplitude(spike_count,1)=temp_sprunghoehe;
            end
        end
    end
end

for cell=1:number_of_cells
        spike_train=find(spike_times_binary(cell,:)==1)';
        estimated_spikes.spike_trains{cell,1}=spike_train;
end
estimated_spikes.raster_plot=spike_times_binary;


if show_figures==1
    figure
    limit=min(number_of_cells, 3);
    for i=1:limit
        eval(sprintf('a%i=subplot(limit,1,%i);',i, i));
        plot(traces_raw(i,:), 'Color', rand(1,3))
        hold on
        peaknumber=1;
        for peaknr=1:size(spike_amplitude_log(i,1).frame,1)
            cellscatter(2,peaknr)=spike_amplitude_log(i,1).frame(peaknr,1);
            cellscatter(1,peaknr)=traces_raw(i,cellscatter(2,peaknr))+1.1.*traces_raw(i,cellscatter(2,peaknr));
            
        end
        scatter(cellscatter(2,:), cellscatter(1,:), 'v')
        title(sprintf('Fluorescence trace with detected peaks of cell %i', i));
        clearvars cellscatter
    end
    linkaxes([a1,a2,a3], 'x')
    
    if number_of_cells>4
        figure
        limit=min(number_of_cells, 6);
        for i=4:limit
            eval(sprintf('a%i=subplot(limit-3,1,%i);',i, i-3));
            plot(traces_raw(i,:), 'Color', rand(1,3))
            hold on
            peaknumber=1;
            for peaknr=1:size(spike_amplitude_log(i,1).frame,1)
                cellscatter(2,peaknr)=spike_amplitude_log(i,1).frame(peaknr,1);
                cellscatter(1,peaknr)=traces_raw(i,cellscatter(2,peaknr))+1.1.*traces_raw(i,cellscatter(2,peaknr));
                
            end
            scatter(cellscatter(2,:), cellscatter(1,:), 'v')
            title(sprintf('Fluorescence trace with detected peaks of cell %i', i));
            clearvars cellscatter
        end
        linkaxes([a4,a5,a6], 'x')
    end
    
    
    figure
    imagesc(spike_times_weighted)
    title('Spiketime raster plot')
    xlabel('image')
    ylabel('cell')
end

fprintf('Finished peak detection\n')
end





