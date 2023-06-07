% i = 6;

normalise_over_stims = 0;%0 for normalising wrt time bins, 1 for normalising wrt stimulation number
all = find(frM.Recording == i & cellfun(@(c) sum(sum(c,1)),frM.Fr_early) >= 10/L);

baseline_avgs = frM.Fr_average_spont(all);
fr_early = frM.Fr_early(all);
fr_middle = frM.Fr_middle(all);
fr_late = frM.Fr_late(all);
% save(sprintf("workspaces/up_down/%s/%i.mat",save_var,i),"frM","L")
% %%
% early_avg = zeros(length(fr_early),1);
% middle_avg = zeros(length(fr_middle),1);
% late_avg = zeros(length(fr_late),1);
% for i = 1:length(fr_early)
%     early_avg(i,:) = mean(mean(base_normalise(baseline_avgs(i),fr_early{i}),normalise_over_stims+1));
%     middle_avg(i,:) = mean(mean(base_normalise(baseline_avgs(i),fr_middle{i}),normalise_over_stims+1));
%     late_avg(i,:) = mean(mean(base_normalise(baseline_avgs(i),fr_late{i}),normalise_over_stims+1));
% end
% %%
% figure
% 
% hold on
% plot([0,1],[0,mean(middle_avg)-mean(early_avg)],"k-","LineWidth",2)
% for i = 1:length(fr_early)
%     if (frM.depth(all(i))<max(frM.depth)/2)
%         plot([0,1],[0,middle_avg(i)./early_avg(i)],"bo-")
%     else
%         plot([0,1],[0,middle_avg(i)-early_avg(i)],"ro-")
%     end
% end
% 
% 

%%
files = dir("workspaces/beforevafter/upregulatedxd/master/*/*.mat");
filed_stims = [];
loop_var =zeros(1,length(files));
for i =1:length(files)
    loop_var(i) = str2num(files(i).folder(111:end));%delete .mat
end
early_avg = cell(max(loop_var),1);
middle_avg = cell(max(loop_var),1);
late_avg = cell(max(loop_var),1);
depths = cell(max(loop_var),1);
max_depth = 0;
unique_depths = [];
for iter_var = loop_var
    iter_var
    if(iter_var ~= 37)
    
    
    load(sprintf("workspaces/up_down/masters/%i.mat",iter_var),"frM","L")
        if any("upregulated" == string(frM.Properties.VariableNames))
            try    
                depths_temp = frM.depth;
                load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",iter_var),"frM","L")
    
                normalise_over_stims = 0;%0 for normalising wrt time bins, 1 for normalising wrt stimulation number
                L = 0.005;
                all = find(frM.Recording == iter_var & (cellfun(@(c) sum(sum(c,1)),frM.Fr_early) >= 1/L | cellfun(@(c) sum(sum(c,1)),frM.Fr_middle) >= 1/L));
                if length(all)>1    
    %                 baseline_stds = frM.baseline_stds(all);
    %                 baseline_avgs = frM.Fr_average_spont(all);
                    fr_early = frM.Fr_early(all);
                    fr_middle = frM.Fr_middle(all);
                    if size(fr_early{1},1) == 1
                        failed_stims = [failed_stims iter_var];
                    end
                    for i = 1:length(fr_early)
                        %%%%%% check if std baseline is 0
    %                     m = mean(reshape([fr_early{i}; fr_middle{i}],1,[]));
    %                     s = std(reshape([fr_early{i}; fr_middle{i}],1,[]));
                        first = 0;%ceil(size(fr_early{i},2)/2);
                        early_avg{iter_var}(i) = mean(reshape(fr_early{i}(:,ceil(size(fr_early{i},2)/2)+1-first:ceil(size(fr_early{i},2))-first),1,[]));%%%%-m/s;%normalise over all stimuli
                        middle_avg{iter_var}(i) = mean(reshape(fr_middle{i}(:,ceil(size(fr_middle{i},2)/2)+1-first:ceil(size(fr_middle{i},2))-first),1,[]));
%                         early_avg{iter_var}(i) = mean(reshape(fr_early{i},1,[]));%%%%-m/s;%normalise over all stimuli
%                         middle_avg{iter_var}(i) = mean(reshape(fr_middle{i},1,[]));
                        depths{iter_var}(i) = depths_temp(all(i));
                        unique_depths = unique([unique_depths;depths{iter_var}(i)]);
                        max_depth = max(max_depth,depths{iter_var}(i));
                    end
                end
            catch
                continue
            end
        end
    end   
        


end

%% Extract activations and dorsal/ventral activations as well as activations for each depth
depth_thr = 1000;
figure
dorsal = cell(length(loop_var),1);%All dorsal neurons of each neuron
ventral = cell(length(loop_var),1);%All ventral neurons of each neuron
activations = cell(length(loop_var),1);%General activations/regulations of all neurons (not dorsal or ventral)
depth_activations = cell(length(unique_depths),1);%All activations ifo their depth
valid_groups = [];

for i = loop_var
    m = mean(middle_avg{i}-early_avg{i});
    s = std(middle_avg{i}-early_avg{i});
    hold on
    %plot([0,1],[0,mean(middle_avg{i}-early_avg{i})],"k-","LineWidth",2)
    activations{i} = (middle_avg{i}-early_avg{i})/s;
    dorsal{i} = (middle_avg{i}(depths{i}<depth_thr)-early_avg{i}(depths{i}<depth_thr))/s;
    ventral{i} = (middle_avg{i}(depths{i}>=depth_thr)-early_avg{i}(depths{i}>=depth_thr))/s;
    if (dorsal{i} ~= 0 & ~isnan(dorsal{i}))
%         plot([0,1],[0,mean(dorsal{i})],"bo-")%look at ventral vs dorsal activation compared to the mean
%         
%         plot([0,1],[0,mean(ventral{i})],"ro-")
        plot((depths{i}(depths{i}<depth_thr)),(dorsal{i}),"b.")
        plot((depths{i}(depths{i}>=depth_thr)),(ventral{i}),"r.")
        valid_groups(end+1) = i;
        for k = 1:length(unique_depths)
            depth_activations{k} = [depth_activations{k} activations{i}(depths{i}==unique_depths(k))];
        end
    end
    
end

legend('','Location', 'northeast')
cols = ["blue","red"];%Source: https://nl.mathworks.com/matlabcentral/answers/334249-manually-enter-legend-details
col_names = ["Dorsal","Ventral"];
for j =1:length(col_names)
    plot([NaN NaN], [NaN NaN], 'Color', cols(j), 'DisplayName', col_names(j))
end
title("Regulation of dorsal & ventral neurons compared to the mean regulation")
ylabel("Mean activation compared to early period")
xlabel("Neuron depth [µm]")
% names = {'Early';'Late'};
% set(gca,'xtick',[0 1],'xticklabel',names)

%% Permutation of the data but in groups 
no_iters = 1000;
means = zeros(no_iters,1);
for iter = 1:no_iters
    perms = randperm(length(valid_groups)*2,length(valid_groups)*2);
    all_activations = cell(length(valid_groups)*2, 1);
    for idx = 1:length(valid_groups)
        all_activations{idx*2-1} = dorsal{valid_groups(idx)};
        all_activations{idx*2} = ventral{valid_groups(idx)};
    end
    dorsals = [];
    ventrals = [];
    for idx = 1:length(all_activations)
        if perms(idx) > length(valid_groups)
            dorsals(end+1:end+length(all_activations{idx})) = all_activations{idx};
        else
            ventrals(end+1:end+length(all_activations{idx})) = all_activations{idx};
        end
    end
    means(iter) = mean(dorsals)-mean(ventrals);
end

%% Compute the unpaired t-test on all neurons
dorsals = [];
ventrals = [];
all_depths = [];
all_neurons = [];
mouse_group = [];
for i = 1:length(dorsal)
    if (dorsal{i} ~= 0 & ~isnan(dorsal{i}))
        dorsals = [dorsals (dorsal{i})];
        ventrals = [ventrals (ventral{i})];
        all_depths = [all_depths depths{i}];
        all_neurons = [all_neurons activations{i}];
        mouse_group = [mouse_group repmat([i],1,length(activations{i}))];
    end

end
[H,P] = ttest2(dorsals,ventrals);

[m,s] = normfit(means);
p_values = tcdf((mean(dorsals)-mean(ventrals)-m)/s,no_iters-1);
%% Create table to fit lme model --> activation is dependent on fixed variable "depth", also on random intercept (and slope?) from mouse group (slope if the variability in a mouse is higher I think) 
tbl = table();
tbl.mouse_group = mouse_group';
tbl.all_neurons = all_neurons';
tbl.all_depths = normalize(all_depths','range');
options = optimset('MaxIter',10000);
% lme = fitlme(tbl, "all_neurons~all_depths+(1|mouse_group)","StartMethod","random","OptimizerOptions",options)
lme2 = fitlme(tbl, "all_neurons~all_depths+(1|mouse_group)+(all_depths-1|mouse_group)")%,"Exclude", [find(residuals(lme) > 2); find(residuals(lme) < -2)])

% compare(lme2,lme)

%% plot activations of middle-early parts of all neurons
figure
for iter = 1:length(dorsal)
    for iter2 = 1:length(dorsal{iter})
    plot([0,1],[0,dorsal{iter}(iter2)],"bo-")%look at ventral vs dorsal activation compared to the mean
    hold on
    end
    
end
for iter = 1:length(ventral)
    for iter2 = 1:length(ventral{iter})
    plot([0,1],[0,ventral{iter}(iter2)],"ro-")%look at ventral vs dorsal activation compared to the mean
    hold on
    end
    
end

 %% change something in all yoke or master tables if needed
% files = dir("workspaces/up_down/yokes/*.mat");
% loop_var =zeros(1,length(files));
% for i =1:length(files)
%     loop_var(i) = str2num(files(i).name(1:end-4));%delete .mat
% end
% for iter_var = [24]
%     iter_var
%     load(sprintf("workspaces/up_down/masters/%i.mat",iter_var),"frM","L")
%     if any("upregulated" == string(frM.Properties.VariableNames))
%         baseline_stds = cellfun(@(c) std(c),frM.Fr_Horridge);
%         frM.baseline_stds = baseline_stds;
%         save(sprintf("workspaces/up_down/masters/%i.mat",iter_var),"frM","L")
%     end
% end
% %  
% %% Create plots of regulations ifo depth
% %source:https://nl.mathworks.com/matlabcentral/answers/180829-shade-area-between-graphs
% depth_std = cellfun(@(c) std(c),depth_activations);
% depth_mean = cellfun(@(c) mean(c),depth_activations);
% curve1 = depth_mean + depth_std;
% curve2 = depth_mean - depth_std;
% x2 = [unique_depths', fliplr(unique_depths')];
% inBetween = [curve1', fliplr(curve2')];
% figure
% fill(x2, inBetween,[1 1 1]*0.9);
% hold on;
% plot(unique_depths, depth_mean, 'k', 'LineWidth', 1);
% 



%% Extract activations and dorsal/ventral activations as well as activations for each depth
depth_thr = 1000;

dorsal_middle = cell(length(loop_var),1);%All dorsal neurons of each neuron
dorsal_early = cell(length(loop_var),1);%All dorsal neurons of each neuron
ventral_middle = cell(length(loop_var),1);%All ventral neurons of each neuron
ventral_early = cell(length(loop_var),1);%All ventral neurons of each neuron
dorsal = cell(length(loop_var),1);%All dorsal neurons of each neuron
ventral = cell(length(loop_var),1);%All ventral neurons of each neuron

activations = cell(length(loop_var),1);%General activations/regulations of all neurons (not dorsal or ventral)
depth_activations = cell(length(unique_depths),1);%All activations ifo their depth
valid_groups = [];
diffs = cell(length(loop_var),1);%All dorsal neurons of each neuron
mean_diffs = [];%All dorsal neurons of each neuron
mean_dorsal = [];
mean_ventral = [];
for i = loop_var
    m = mean(middle_avg{i}-early_avg{i});
    s = std(middle_avg{i}-early_avg{i});
    %plot([0,1],[0,mean(middle_avg{i}-early_avg{i})],"k-","LineWidth",2)
    activations{i} = (middle_avg{i}-early_avg{i})/s;
    dorsal_middle{i} = middle_avg{i}(depths{i}<depth_thr);
    dorsal_early{i} = early_avg{i}(depths{i}<depth_thr);
    ventral_middle{i} = middle_avg{i}(depths{i}>=depth_thr);
    ventral_early{i} = early_avg{i}(depths{i}>=depth_thr);
    dorsal{i} = (middle_avg{i}(depths{i}<depth_thr)-early_avg{i}(depths{i}<depth_thr));
    ventral{i} = (middle_avg{i}(depths{i}>=depth_thr)-early_avg{i}(depths{i}>=depth_thr));
    if ~(isempty(dorsal{i}))
        mean_dorsal(end+1) = mean(dorsal{i});
    end
    if ~(isempty(ventral{i}))
        mean_ventral(end+1) = mean(ventral{i});
    end
    if ~(isempty(dorsal{i}) || isempty(ventral{i}))
        mean_diffs(end+1) = mean(dorsal{i})-mean(ventral{i});
        
    end
    
end
%% bootstrap resampling
reps = 10000;
means_master = zeros(1,reps);
for iter = 1:reps
    idxs = randi([1 length(mean_diffs)],length(mean_diffs),1);
    means_master(iter) = mean(mean_diffs(idxs),"omitnan");
end

figure
histogram(means_master)
signrank(mean_diffs)
