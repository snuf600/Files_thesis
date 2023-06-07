clear
close all

load("memory_cells\height_spike_rate\7.mat")
L = 1
i = recording_number;
plotting = 1;
normalise_over_stims = 0;%0 for normalising wrt time bins, 1 for normalising wrt stimulation number
all = find(frM.Recording == i & cellfun(@(c1,c2) sum(sum(c1(~isnan(c1)),1)),frM.Fr_horridge,frM.heights) >= 100/L);

heights = frM.heights(all);
fr_horridge = frM.Fr_horridge(all);


meanning = 0;
if meanning
    for k = 1:length(all)
        heights{k} = (heights{k}(1:2:end-1)+heights{k}(2:2:end))/2;
        fr_horridge{k} = (fr_horridge{k}(1:2:end-1)+fr_horridge{k}(2:2:end))/2;
    end
end



if plotting
    for k = 1:length(all)

%         heights{k}(1:70) = nan;
        figure
        col = 0:1:(length(heights{k})-1);
        scatter(fr_horridge{k},heights{k},7.5,col,"filled")
        xlabel("Paw Height [mm]")
        ylabel("Firing rate [Hz]")
        title("Relation between paw height and firing rate")
        c = colorbar;
        c.Label.String = "Time [s]";
        caxis([find(~isnan(heights{k})==1,1,"first") find(~isnan(heights{k})==1,1,"last")])
        %Calculate mean/height bin
%         edges = (0:0.1:1)*(max(heights{k})-min(heights{k}))+min(heights{k});
%         means = zeros(length(edges)-1,1);
%         mean_heights = zeros(length(edges)-1,1);
%         for i = 1:length(edges)-1
%             means(i) = mean(fr_horridge{k}(heights{k} < edges(i+1) & heights{k} > edges(i)));
%             mean_heights(i) = sum(edges(i:i+1))/2;
%         end
%         mean_heights = mean_heights(~isnan(means));
%         means = means(~isnan(means));
%         interp_heights = min(mean_heights):0.001:max(mean_heights);
%         means_interpolated = spline(mean_heights,means,interp_heights');
%         hold on
%         yyaxis right
% %         plot(interp_heights,means_interpolated,"r")
%         plot(mean_heights,means,"r")
%         ylim([0 max(means)])
    end
end

figure
plot(heights{1})
