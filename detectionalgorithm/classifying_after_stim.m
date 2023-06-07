% %% Generate valid points data by looking at no stimulation data
% % load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220\events.mat')
% % load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220831\t.mat')
% % load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220831\ZmasterD1P2.mat')
% zData = ZMaster(:,1);
% onsets = events.onsets(events.onsets<1600 & events.onsets>900)-900;
% fs = 1/L;
% moving_window = round(fs*2);%2s without stimulations
% %Compute points that don't have spikes or random weird noise for window of 10s
% valid_points = ones(length(zData),1);%Valid where no spike within window
% for i = 1:length(onsets)
%     %First find allocated datapoint to timestamp
%     [~,idx] = min(abs(t-onsets(i)));
%     %Use this index to set window around it to 0
%     valid_points(max(idx-moving_window,1:min(length(valid_points),idx+moving_window))) = 0;
% end
% % Filter out real noisy parts
% moving_std = movstd(zData,100);
% moving_std_thr = 0.20;
% thresholded_std = moving_std > moving_std_thr;
% valid_points_std = ones(length(zData),1);
% for i = 1:length(zData)
%     if i<moving_window
%         valid_points_std(i) = sum(ismember(onsets,1:i+moving_window));
%     elseif i>length(zData)-moving_window
%         valid_points_std(i) = sum(thresholded_std(i-moving_window+1:end));
%     else
%         valid_points_std(i) = sum(thresholded_std(i-moving_window+1:i+moving_window));
%     end
% end
% valid_points_std = valid_points_std == 0;%Where sum == 0, no spikes are nearby
% valid_points = valid_points_std.*valid_points;
% 
% 
% %Take out stimulation contaminated data
% valids = unique(round(t(valid_points == 1)/L));
% valids = [valids [max(valids)+1:600/L]];%After last datapoint all points are valid
% % baseline_stds = zeros(length(all),1);
% % baseline_avgs = zeros(length(all),1);
% % for i = 1:length(all)
% %     baseline_stds(i) = std(fr_horridge_matrix(i,valids));
% %     baseline_avgs(i) = mean(fr_horridge_matrix(i,valids));
% % end

%% Load data
i = recording_number;
plotting = 1;
normalise_over_stims = 0;%0 for normalising wrt time bins, 1 for normalising wrt stimulation number
all = find(frM.Recording == i & cellfun(@(c) sum(sum(c,1)),frM.Fr_early)+cellfun(@(c) sum(sum(c,1)),frM.Fr_middle) >= 100/L);

fr_early = frM.Fr_early(all);
fr_middle = frM.Fr_middle(all);
fr_horridge = frM.Fr_Horridge(all);
fr_horridge_matrix = zeros(length(fr_horridge),length(fr_horridge{1}));
early = zeros(length(fr_early),size(fr_early{1},1),size(fr_early{1},2));
middle = zeros(length(fr_middle),size(fr_middle{1},1),size(fr_middle{1},2));
%Generate matrix of whole horridge period
for i = 1:length(all)
    fr_horridge_matrix(i,:) = fr_horridge{i}*L; 
    early(i,:,:) = fr_early{i}*L;
    middle(i,:,:) = fr_middle{i}*L;
end
activations = [early middle];

% Generate avgs over stimulation
avgs = zeros(length(fr_early),size(fr_early{1},(-1)*normalise_over_stims+2));

% To store the z-scores of all neurons compared to random intervals
p_values = zeros(length(fr_early),size(fr_early{1},(-1)*normalise_over_stims+2));

for i = 1:length(fr_early)
    avgs(i,:) = normalize((mean(reshape(early(i,:,:),size(early,2),[]),normalise_over_stims+1)*size(early,2) ...
        +mean(reshape(middle(i,:,:),size(middle,2),[]),normalise_over_stims+1)*size(middle,2)) ...
        /(size(early,2)+size(middle,2)));
%     avgs(i,:) = normalize(avgs(i,:),"zscore","robust");
end
class = cell(length(all),1);
%% Create a sliding window to check where a change in time occurs then use this period in which the change occurs to check whether or not there is also a change through stimulation number
width = 6;
window = [ones(width/2,1); zeros(width/2,1)]-1/2;
p_values = zeros(length(all),1);
start_end = zeros(length(all),2);
thr = 10; %avg of 3 z-scores higher/lower
relevant_neurons = [];
for j = 1:length(all)
%     if j == 2
%         aaaaaa = 1
%     end
    activation = avgs(j,floor(length(avgs(j,:))/2)+1:end);%1x40
    
%     midpoint = ceil(length(activation)/2);
%     convoluted_activation = conv2(activation',window,"valid");%1x35, from middle on (18th value until end to check max difference between before and after)
    
    step_size = 2;
    total_sum = sum(activation,"all");
    max_range_sum = 0;
    for range = 1:length(activation)-step_size
        if sum(activation(range:range+step_size)) > max_range_sum
            max_range = [range range+step_size];
            max_range_sum = sum(activation(range:range+step_size));
        end
    end
    
    step_size = 2;
    total_sum = sum(activation,"all");
    min_range_sum = max_range_sum;
    for range = 1:length(activation)-step_size
        if sum(activation(range:range+step_size)) < min_range_sum
            min_range = [range range+step_size];
            min_range_sum = sum(activation(range:range+step_size));
        end
    end
    if abs(min_range_sum) < abs(max_range_sum)
        val = abs(max_range_sum);
    else
        val = (min_range_sum);
    end
    
    if val > 5%on avg 3std below zero
        relevant_neurons(end+1) = j;
        start_end(j,:) = max_range+floor(length(avgs(j,:))/2);
    elseif val < -5
        class{j} = "B";
    else
        class{j} = "None";
    end
    
%     start_after = ceil(length(convoluted_activation)/2);
%     before_stim = convoluted_activation(1:start_after-1);
%     after_stim = convoluted_activation(start_after:end);
%     
%     [min_value,idx1] = min(after_stim);%activity goes down from this timebin on
%     [max_value,idx2] = max(after_stim);%activity goes up from this timebin on
%     if (max(abs(min_value),abs(max_value)) > thr && idx2<idx1)
%         start_period = min(idx1,idx2);
%         end_period = min(max(idx1,idx2),start_period+4);%max 5 timebins
%         start_end(j,:) = [start_period+midpoint end_period+midpoint];
%         relevant_neurons(end+1) = j;
%     end
end

%% Check whether there is a change through stimulation number
nb_iters = 1000;
p_value = zeros(length(relevant_neurons),1);
p2 = zeros(length(relevant_neurons),1);
stim_number = zeros(length(relevant_neurons),1);
min_window_size = ceil(0.15*size(activations,2));%Values in beginning are ignored due to high noise sometimes
% window_size = ceil(0.1*size(activations,2));%Determines the bins around the optimal threshold that are used to compute the difference
plotting = 0;
if master
    mkdir(sprintf("./workspaces/beforevafter/%i/master/%i",upup,recording_number))
else
    mkdir(sprintf("./workspaces/beforevafter/%i/yoke/%i",upup,recording_number))
end
for iter = relevant_neurons
    
    disp(sprintf("Master %i, neuron %i",recording_number,iter))
    data = reshape(activations(iter,:,start_end(iter,1):start_end(iter,2)),size(activations,2),[]);%Select correct period from a certain neuron
    data_mean = mean(data,2);
    [idx,fun] = find_best_change(data_mean,min_window_size,upup);
%     idx = ceil(length(data_mean)/2);
    window_size = min(idx-1,length(data_mean)-idx);%max possible window size
    max_diff = (mean(data_mean(idx:idx+window_size-1))-mean(data_mean(idx-window_size:idx-1)));
    dist = zeros(nb_iters,1);%to store the random permutated values
    for i = 1:nb_iters%compute same thing for permuted versions of the same array
        arr_perm = zeros(size(data));
        for k = 1:size(arr_perm,2)
            arr_perm(:,k) = data(randperm(size(data,1)),k);
        end
        arr_perm = mean(arr_perm,2);
        [~,idx_perm] = max(arrayfun(@(i) abs(mean(arr_perm(i:i+window_size-1))-mean(arr_perm(i-window_size:i-1))),window_size+1:length(arr_perm)-window_size+1));%compute best stim number threshold
        idx_perm = idx_perm + window_size;
%             max_diff_perm = (mean(arr_perm(idx:end))-mean(arr_perm(1:idx-1)));
        max_diff_perm = abs(mean(arr_perm(idx_perm:idx_perm+window_size-1))-mean(arr_perm(idx_perm-window_size:idx_perm-1)));
        dist(i) = max_diff_perm;
    end
    p_value(relevant_neurons == iter) = (sum(dist>=abs(max_diff))+1)/(nb_iters+1);%check what fraction of permuted samples are at least as unlikely as the one recorde
    stim_number(relevant_neurons == iter) = idx;
    p2(relevant_neurons == iter) = c_test(data(1:idx,:),data(idx+1:end,:));%if neuron upregulated --> close to 1, if downregulated close to 0. Compares two distributions separated by idx
    if p_value(relevant_neurons == iter) <= 0.005
        if (p2(relevant_neurons == iter) <= 0.001 && ~upup)
            class{iter} = "C";
        elseif (1-p2(relevant_neurons == iter) <= 0.001 && upup)
            class{iter} = "D";
        else
            class{iter} = "A";
        end
    else
        class{iter} = "A";
    end

    if plotting    
        fig = figure("visible","off");
    %     hist(dist)
        imagesc(reshape(activations(iter,:,:),size(activations,2),[]));
        title(sprintf("Neuron %i, class %s, p1 value = %.3f,p2 value = %.3f",iter,class{iter},p_value(relevant_neurons == iter),p2(relevant_neurons == iter)))
        
            yline(idx,"r","LineWidth",4)
        
%         subplot(2,1,2)
% 
%         plot(1+min_window_size:length(fun)+min_window_size,fun);
%         hold on
%         plot(data_mean,"r")
    
    if master
        saveas(fig,sprintf("./workspaces/beforevafter/%i/master/%i/%i.png",upup,recording_number,iter));
    else
        saveas(fig,sprintf("./workspaces/beforevafter/%i/yoke/%i/%i.png",upup,recording_number,iter));
    end
    end
    close all
    
end
if master
    save(sprintf("./workspaces/beforevafter/%i/master/%i/data.mat",upup,recording_number),"events","frM","L","p_value","stim_number","p2","relevant_neurons","class");
else
    save(sprintf("./workspaces/beforevafter/%i/yoke/%i/data.mat",upup,recording_number),"events","frM","L","p_value","stim_number","p2","relevant_neurons","class");
end
%%
for i = 1:length(all)'
    if ~any((ismember(relevant_neurons,i)))
        if plotting    
            fig = figure("visible","off");
        %     hist(dist)
            imagesc(reshape(activations(i,:,:),size(activations,2),[]));
            title(sprintf("Neuron %i, class %s",i,class{i}))
            
            if master
                saveas(fig,sprintf("./workspaces/beforevafter/%i/master/%i/%i.jpg",upup,recording_number,i));
            else
                saveas(fig,sprintf("./workspaces/beforevafter/%i/yoke/%i/%i.jpg",upup,recording_number,i));
            end
            close all
        end
    end

end





%% Find the fraction of significant downregulated neurons
aaa = 0;
if aaa
    clear
    recording_database_updated_20200805_JCfreshStart
    
    sort_horridge = sort_masters_horridge;
    sort_before = sort_masters_before;
        
    
    significant_down_regulation_master = 0;
    total_neurons_master = 0;
    
    for i = 1:length(sort_horridge)%; 5 is 23/09
        
        try
            i
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
            relevant_neurons((p_value<=0.001 & p2 < 0.001));
            significant_down_regulation_master = significant_down_regulation_master + sum(p_value<=0.001 & 1-p2 < 0.001);
            total_neurons_master = total_neurons_master + length(relevant_neurons);%length(frM.Fr_early);%%%%%%%
    
        catch
            continue
        end
    end
    
    
    sort_horridge = sort_yokes_horridge;
    sort_before = sort_yokes_before;
    
    significant_down_regulation_yoke = 0;
    total_neurons_yoke = 0;
    
    for i = 1:length(sort_horridge)%; 5 is 23/09
        
        try
            i
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
            (p_value<=0.001 & p2 < 0.001);
            significant_down_regulation_yoke = significant_down_regulation_yoke + sum(p_value<=0.001 & 1-p2 < 0.001);
            total_neurons_yoke = total_neurons_yoke + length(relevant_neurons);%length(frM.Fr_early);%%%%%%%
    
    
        catch
            continue
        end
    end
    p1 = significant_down_regulation_master/total_neurons_master; p2 = significant_down_regulation_yoke/total_neurons_yoke; n1 = total_neurons_master;n2=total_neurons_yoke;p=(significant_down_regulation_master+significant_down_regulation_yoke)/(n1+n2);
z = (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
end
%% Find the fraction of significant downregulated neurons
aaa = 0;
if aaa
    clear
    recording_database_updated_20200805_JCfreshStart
    
    sort_horridge = sort_masters_horridge;
    sort_before = sort_masters_before;
        
    
    significant_down_regulation_master = 0;
    total_neurons_master = 0;
    probs_master = [];
    success_total_master = [];
    masters = 0;
    for i = 1:length(sort_horridge)%; 5 is 23/09
        
        try
            
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
            relevant_neurons((p_value<=0.001 & p2 < 0.001));
            significant_down_regulation_master = significant_down_regulation_master + sum(p_value<=0.001 & 1-p2 < 0.001);
            total_neurons_master = total_neurons_master + length(frM.Fr_early);%%%%%%%
            probs_master(end+1) = sum(p_value<=0.001 & 1-p2 < 0.001)/length(relevant_neurons);
            success_total_master(end+1,:) = [sum(p_value<=0.001 & 1-p2 < 0.001) length(frM.Fr_early)];
%             if sum(p_value<=0.001 & 1-p2 < 0.001)>0
%                 masters = masters+1;
%             end
%             load(sprintf("workspaces/beforevafter/1/master/%i/data.mat",i))
%             probs_master(end+1) = sum(p_value<=0.001 & 1-p2 < 0.001)/length(relevant_neurons);
% 
%             success_total_master(end,:) = success_total_master(end,:)+ [sum(p_value<=0.001 & 1-p2 < 0.001) 0];
            
        catch
            continue
        end
    end
    
    
    sort_horridge = sort_yokes_horridge;
    sort_before = sort_yokes_before;
    
    significant_down_regulation_yoke = 0;
    total_neurons_yoke = 0;
    
    probs_yoke = [];
    success_total_yoke = [];
    yokes = 0;
    for i = 1:length(sort_horridge)%; 5 is 23/09
        
        try
            
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
            (p_value<=0.001 & p2 < 0.001);
            significant_down_regulation_yoke = significant_down_regulation_yoke + sum(p_value<=0.001 & 1-p2 < 0.001);
            total_neurons_yoke = total_neurons_yoke + length(frM.Fr_early);%%%%%%%
            probs_yoke(end+1) = sum(p_value<=0.001 & 1-p2 < 0.001)/length(relevant_neurons);

            success_total_yoke(end+1,:) = [sum(p_value<=0.001 & 1-p2 < 0.001) length(frM.Fr_early)];
%             if sum(p_value<=0.001 & 1-p2 < 0.001)>0
%                 yokes = yokes+1;
%             end
%             load(sprintf("workspaces/beforevafter/1/yoke/%i/data.mat",i))
%             probs_yoke(end+1) = sum(p_value<=0.001 & 1-p2 < 0.001)/length(relevant_neurons);
% 
%             success_total_yoke(end,:) = success_total_yoke(end,:) + [sum(p_value<=0.001 & 1-p2 < 0.001) 0];

        catch
            continue
        end
        
    end
    p1 = significant_down_regulation_master/total_neurons_master; p2 = significant_down_regulation_yoke/total_neurons_yoke; n1 = total_neurons_master;n2=total_neurons_yoke;p=(significant_down_regulation_master+significant_down_regulation_yoke)/(n1+n2);
    z = (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
end
% %% Permutation of the data but in groups 
% success_total = [success_total_master ; success_total_yoke];
% shuffling_vector = [ones(1,length(success_total_master)) zeros(1,length(success_total_yoke))];
% no_iters = 10000;
% diffs = zeros(no_iters,1);
% for iter = 1:no_iters
%     shuffled_vector = shuffling_vector(randperm(length(shuffling_vector)));
%     mean_p_master = mean(success_total(shuffled_vector == 1,1)./success_total(shuffled_vector == 1,2),"omitnan");
%     mean_p_yoke = mean(success_total(shuffled_vector == 0,1)./success_total(shuffled_vector == 0,2),"omitnan");
%     diffs(iter) = mean_p_master-mean_p_yoke;
% end
% real_mean_p_master = mean(success_total_master(:,1)./success_total_master(:,2),"omitnan");
% real_mean_p_yoke = mean(success_total_yoke(:,1)./success_total_yoke(:,2),"omitnan");
% real_diff = real_mean_p_master-real_mean_p_yoke;
% %% bootstrap resampling
% reps = 10000;
% means_master = zeros(1,reps);
% for iter = 1:reps
%     idxs = randi([1 length(probs_master)],length(probs_master),1);
%     means_master(iter) = mean(probs_master(idxs),"omitnan");
% end
% reps = 10000;
% means_yoke = zeros(1,reps);
% for iter = 1:reps
%     idxs = randi([1 length(probs_yoke)],length(probs_master),1);
%     means_yoke(iter) = mean(probs_yoke(idxs),"omitnan");
% end
% figure
% histogram(means_master-means_yoke)
% [h,p] = ttest2(means_master,means_yoke)
%%

 function [idx,fun] = find_best_change(data_mean,ignore,up_reg)
    data_mean_norm = normalize(data_mean);
    optimise_function = arrayfun(@(i) (mean(data_mean(i:end))-mean(data_mean(1:i-1))),1+ignore:length(data_mean)-ignore);%compute best stim number threshold
    ema_filtered = 1+normalize(abs(arrayfun(@(i) sum(data_mean_norm(1:i)),1:length(data_mean_norm))),"range");
    %     ema_filtered = zeros(1,length(data_mean));
%     for iter = 2:length(data_mean)
%         if data_mean(iter)>0
%             ema_filtered(iter) = 0.95*ema_filtered(iter-1)+0.05*0.2;
%         else
%             ema_filtered(iter) = 0.95*ema_filtered(iter-1);
%         end
%     end
    fun = [optimise_function];
    optimise_function = optimise_function.*(ema_filtered(1+ignore:length(data_mean_norm)-ignore));
    if up_reg   
    [~,idx] =  max(optimise_function);

    else
    [~,idx] =  min(optimise_function);
    end
    fun = [fun;ema_filtered(1+ignore:length(data_mean_norm)-ignore)];
    idx = idx+ignore;
   
    
end       








% %% Create a distribution of the mean and std of every neuron for the amount of stimulations
% reps = 1000;
% means = zeros(length(all),reps);
% m = zeros(1,length(all));
% s = zeros(1,length(all));
% for j = 1:length(all)
% %     all(j)
%     for iter = 1:reps
%         idxs = randperm(length(valids),size(activations,2));%create a matrix containing which random sample points to use
%         random_snippets = zeros(1,size(activations,2));
%         random_snippets = fr_horridge{j}(valids(idxs))*L;
%         means(j,iter) = mean(random_snippets)*size(activations,2);%poisson distributed process
%     end
%     lambda = poissfit(means(j,:));
%     plotting = 1;
%     if plotting
%         x = round(unique(means(j,:)));
%         figure;
%         hold on
%         histogram(means(j,:),"BinEdges",min(means(j,:))-0.5:1:max(means(j,:))+0.5,"Normalization","pdf");
%         plot(x,poisspdf(x,lambda))
%         legend("measured means","fitted poisson model")
%     end
%     p_values(j,:) = poisscdf(avgs(j,:),lambda);
% end
% 
% %%
% figure
% imagesc(p_values)
% colorbar
% 
% %%
% p_values(j) = tcdf((means_true(j)-m)/s,reps-1);
% significant_increase = find(p_values>0.975);%two tailed t-test
% p_values(p_values>0.5) = 1-p_values(p_values>0.5);
% p_values = 2*p_values;%2 tailed ttest
% significant_neurons = find(p_values<0.05);
% 
% 
% %% Plot data
% if plotting
%     for k = significant_increase'%1:length(all)
%         figure
%         imagesc(base_normalise(mean(fr_horridge{k}),std(fr_horridge{k}),fr_lifting{k}))
%         xlabel(sprintf("Timebins of %ims bins",L*1000))
%         ylabel("Stimulation number")
%         title(sprintf("Activation of neuron %i in %ims intervals before a stimulation",k,T*1000))
%         colorbar
%     end
% %     figure
% %     imagesc(lifting(significant_increase,:))
% %     xlabel(sprintf("Timebins of %ims bins",L*1000))
% %     ylabel("Neuron number")
% %     title(sprintf("Mean activation of neurons in %ims intervals before a stimulation during lifting period",T*1000))
%     figure
%     imagesc(lifting_avg(significant_increase,:))
%     xlabel(sprintf("Timebins of %ims bins",L*1000))
%     ylabel("Neuron number")
%     title(sprintf("Mean activation of neurons in %ims intervals before a stimulation during lifting period",T*1000))
%     
% %     figure
% %     imagesc(middle)
% %     xlabel(sprintf("Timebins of %ims bins",L*1000))
% %     ylabel("Neuron number")
% %     title(sprintf("Mean activation of neurons in %ims intervals before a stimulation during middle period",T*1000))
% %     
% %%
%     figure
%     subplot(2,1,1)
%     imagesc(fr_horridge_matrix)
%     xlabel(sprintf("Timebins of %ims bins",L*1000))
%     ylabel("Neuron number")
%     title(sprintf("Mean activation of neurons in %ims intervals before a stimulation during lifting period",T*1000))
%     subplot(2,1,2)
%     plot(timestamps,zeros(1,length(timestamps)),"r*")
%     xlim([tbin(2,1) tbin(4,2)])
% %     figure
% %     imagesc(middle_avg)
% %     xlabel(sprintf("Timebins of %ims bins",L*1000))
% %     ylabel("Neuron number")
% %     title(sprintf("Mean activation of neurons in %ims intervals before a stimulation during middle period",T*1000))
% end
% % TO DO: normalise with baseline, 
% 
% %Wat dit zou kunnen zijn: de neuron leert sneller en sneller te reageren op
% %de stimuli piempampoem
% 
% % figure
% % subplot(2,1,1)
% % imagesc([lifting_avg middle_avg])
% % subplot(2,1,2)
% % yyaxis left
% % plot(spike_spacing)
% % xlim([1 102])
% % hold on
% % 
% % segments = zeros(length(spike_spacing),1);
% % mm = mean([lifting middle],1);
% % for i = 1:length(spike_spacing)
% %     segments(i) = mean(mm((i-1)*round(T/L)+1:i*round(T/L)));
% % end
% %yyaxis right
% %plot(segments)
% 
% 
% % [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = ica_assembly_detection(lifting',1);
% % assemblies
