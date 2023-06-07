% clear

recording_database_updated_20200805_JCfreshStart



% i = number of the recording
i = 5;
master = 1;%masters or yokes
if master
    sort_horridge = sort_masters_horridge;
    sort_before = sort_masters_before;
else
    sort_horridge = sort_yokes_horridge;
    sort_before = sort_yokes_before;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Creation of frY and frM table for analysis
%
% 29/04/22
% Simon Lavaud
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% - - - Help on how to use it:
%
% You will probably already have a frM frY table from preivous recording
% The best is to load your previous table, rename it as frM_old frY_old.
% Then you can run for a specific recording the frM frY algorithm.
% It will return a frM frY only for one recording.
% You can copy and past this new frM frY in the old version.
% Then delete the new frM and frY, and rename the frM_old and frY_old
% as frM and frY. You can now save it in the right folder.



% Select 1 for bad masters/yokes and 0 for good ones
% BAD = 0;
% Select the database of interest
% recording_database_updated_20200805_JCfreshStart;

% Initialize the variable
nch = 385;
Fs = 30000;
analysis_type = 1;


% Column of the table for analysis
Recording = [];
Id_spont = [];
Fr_spont = [];
Id_horr = [];
Fr_early = [];
Fr_middle = [];
Fr_late = [];
Fr_waiting = [];
class = [];

% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - frM
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

clearvars frM

% Initialize inner variable
cont = 0;

% Iteration on the recording
% If only one recording to add : i = recording number
% If all the recording: [7 8 10 11 12 13 14 15 16 17 18 19 20 21 22]
%%
for i = 7:7%5:7%1:length(sort_horridge); 5 is 23/09
%     sort_horridge = "W:\Charlotte\Spinal_Cord\Sorted\20221123\horridge"
    clearvars -except i sort_horridge master sort_before spk
    recording_number = i;
    switch i
        case 2
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220819\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220819\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220819\ZmasterD1P1.mat')
        case 3
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220831\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220831\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220831\ZmasterD1P2.mat')
        case 4
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220921\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220921\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220921\ZmasterD1P3.mat')
        
        case 5
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220923\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220923\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\20220923\ZmasterD1P4.mat')
        case 6
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\16_11_2022\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\16_11_2022\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\16_11_2022\ZmasterD1P4.mat')
        case 7
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\23_11_2022\events.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\23_11_2022\t.mat')
        load('D:\From c disk\Documenten\Arnaud\School\5e jaar\NERF\files train\memory_cells\AnalyseDay1 - MASTER\AnalyseDay1 - MASTER\23_11_2022\ZmasterD1P2.mat')
    end
    heights = find_heights(ZMaster(:,1),t,events);
    t_heights = t+900;
    recording_database_updated_20200805_JCfreshStart


%for i = [2:4 7:9]
    
    % Initialize the variable
    nch = 385;
    Fs = 30000;
    analysis_type = 1;

    % Column of the table for analysis
    Recording = [];
    Id_spont = [];
    Fr_spont = [];
    Id_horr = [];
    Fr_early = [];
    Fr_middle = [];
    Fr_late = [];
    Fr_waiting = [];
    class = [];

    recording_number = i;

    % Display the recording that is analysed
    if master
        disp(sprintf('Master %d',i))
    else
        disp(sprintf('Yoke %d',i))
    end
    
    spk = loadKSdir(sort_horridge{i});
    clusters = spk.cids(spk.cgs == 2);
    
    %crete time binning
    
    load(fullfile(sort_horridge{i},'events.mat'),'events');
%     save(sprintf("workspaces/%i.mat",i))

%     load(sprintf("workspaces/%i.mat",i))
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
    end
    
    %Load first and last event for Eearly and Late
    %For it: look at the event file.
    if i==1 %25/07/22 // Ptf1a
        lim_start=1;
        lim_event=324;

    elseif i==2 %19/08/22 // Ptf1a
        lim_start=1;
        lim_event=336;
    
    elseif i==3 %19/08/22 // Eng (waiting: 3676-7264)
        lim_start=1;
        lim_event=749;
        waiting_start = 3676;
        waiting_stop = 7264;
    elseif i==4 %21/09/22 // Ptf1a (waiting: 3709-7382)
        lim_start=2;
        lim_event=1070;
        waiting_start = 3709;
        waiting_stop = 7382;
    elseif i==5 %23/09/22 // Ptf1a 
        lim_start=1;
        lim_event=128;
        waiting_start = 3921;
        waiting_stop = 7924;
    elseif i==6 %11/10/22 // Ptf1a 
        lim_start=1;
        lim_event=103;
        waiting_start = 3940;
        waiting_stop = 7563;
    elseif i==7 %16/11/22 // Ptf1a 
        lim_start=3;
        lim_event=419;
        waiting_start = 3025;
        waiting_stop = 6902;
    elseif i==8 %24/11/22 // ENG 
        lim_start=1;
        lim_event=164;
        waiting_start = 3075;
        waiting_stop = 7212;
    elseif i==9 %28/11/22 // ENG
        lim_start=1;
        lim_event=225;
        waiting_start = 3207;
        waiting_stop = 7213;
        
        
        
        
        % Initialization parameter in case of problem
    else
        lim_start=2;
        assert(events.onsets(lim_start)>=900)
        while events(evtid).onsets(lim_start+1)>events(evtid).onsets(lim_start)+1
            lim_start = lim_start+1;
            if lim_start == 100
                break
            end
        end
        lim_event=sum(events(evtid).onsets(2) <= events(evtid).onsets & events(evtid).onsets <= events(evtid).onsets(2)+600)+1;
    end
    
    
    if lim_start<1
        lim_start=2;
    end
    
    % Computation of the timeline Early, mid, late
    % Mid not use anymore
    clearvars dstim tdist
    a=1;
    %Frequency
    for j = (lim_start+1):lim_event
        tdist(j) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
        %a=a+1;
    end
    
    dstim = smoothdata(tdist);
    idmax = find(dstim>1,1);
    if isempty(idmax)
        idmax=lim_start+10;
    end
    
    if i==14 | i==16
        filter_non_onset=find(events.onsets>800 & events.onsets<1000,1);
    else
        filter_non_onset=1;
    end
    
    % Tbin
    % 0 to 600 - spont
    % 600 to time early (early phase)
    % time early to time mid (mid phase)
    % time mid to end horridge (late phase)
    % end horridge - end+600 (after phase)
    % waiting phase (50 minutes after horridge): enter new values foreach
    % recording
    
    if ~isempty(idmax)
        tbin = [events(evtid).onsets(lim_start)-600, events(evtid).onsets(lim_start)-1;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(idmax);... %Early
            events(evtid).onsets(idmax), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %late

%             waiting_start, waiting_stop]; %late
        
        % events(evtid).onsets(lim_start)+600, events(evtid).onsets(lim_start)+600+2700; %Between
    else
        tbin = [events(evtid).onsets(1)-600, events(evtid).onsets(1)-1;...
            events(evtid).onsets(1), 60;...
            60, events(evtid).onsets(end);...
            events(evtid).onsets(end), events(evtid).onsets(1)+600];
%             waiting_start,waiting_stop]; 
    end

    assert(tbin(3,2)<tbin(4,2))
    disp(tbin)
    
    %%
    % Create columm :
    % Recording number
    % Id horridge / Id spont / Id after (same on long recording)
    L = 1; %change to increase/decrease time bins
    step = 0.1;
    cont = 0;
    spike_spacing = [];
    long_pause = [];
    

    for j = 1:length(clusters)
        cont = cont+1;
%         if master
%             disp(sprintf('master %d, cluster long %d',i,j))
%         else
%             disp(sprintf('yoke %d, cluster long %d',i,j))
%         end
        Recording(cont) = i;
        Id_horr(cont) = clusters(j);
        Id_spont(cont) = clusters(j);
        Id_after(cont) = clusters(j);
        
        % Create the FR column for each phases d epending on Tbin
        for k = 2
            edges = tbin(k,1):L:tbin(k+2,2);
            mean_heights = zeros(length(edges)-1,1);
            
            for time_var = 1:length(edges)-1
                valid_times = (t_heights>edges(time_var)) & (t_heights<edges(time_var+1));
                mean_heights(time_var) = mean(heights(valid_times),"omitnan");
            end
            
            edges = tbin(k,1):L:tbin(k+2,2);
            frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
            Fr_horridge{cont} = frmov;
            neuron_height{cont} = mean_heights;

            
        end
        k = 1;
        edges = tbin(k,1):L:tbin(k,2);
        frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
        Fr_spont{cont} = frmov;
%         
%         k=5;
%         win = [1500 2100];
%         frmov = [];
%         edges = win(1):L:win(2);
%         frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
%         while win(2)<=2100
%             stj = spk.st(spk.clu == clusters(j) & spk.st>win(1) & spk.st<=win(2));
%             frmov = [frmov numel(stj)/L];
%             win = win+L;
%         end
%         Fr_after{cont} = frmov;
        
%         Fr_average_after(cont)=mean(frmov);
        
        Fr_average_spont(cont)=mean(Fr_spont{cont});
        
%         Fr_average_waiting(cont) = mean(Fr_waiting{cont});

        Tbin_tot{cont}=tbin;
    end



%%
% Create the frM table


frM = table();
frM.Recording = Recording';
frM.Tbin = Tbin_tot';
frM.Id_spont = Id_spont';
frM.Fr_spont = Fr_spont';
frM.Fr_average_spont = Fr_average_spont';
frM.Id_horr = Id_horr';
frM.Fr_horridge = Fr_horridge';
% frM.Fr_pred = Fr_early_stim';
frM.heights = neuron_height';
% frM.Fr_middle_stim = Fr_middle_stim';
% frM.Fr_middle_no_stim = Fr_middle_no_stim';
% frM.Fr_after = Fr_after';
% frM.Fr_average_after = Fr_average_after';
% frM.Fr_waiting = Fr_waiting';
% frM.Fr_average_waiting = Fr_average_waiting';

%%

%firing rates means for 5 bins (10min) during wainting
% Fr_changes_test = zeros(length(Fr_spont),5);
% 
% time_bin = 600*5; %time window for data extraction
% 
% for n = 1:length(Fr_spont)
% 
%     for c = 1:5
%     
% 	    Fr_changes_test(n,c) = mean(Fr_waiting{1,n}(1,((c-1)*time_bin)+1:(time_bin +(c-1)*(time_bin-1))+1));
%     
%     end
% end
% 
% Fr_average_waiting_1 = [];
% Fr_average_waiting_1 = Fr_changes_test(:,1);
% Fr_average_waiting_2 = [];
% Fr_average_waiting_2 = Fr_changes_test(:,2);
% Fr_average_waiting_3 = [];
% Fr_average_waiting_3 = Fr_changes_test(:,3);
% Fr_average_waiting_4 = [];
% Fr_average_waiting_4 = Fr_changes_test(:,4);
% Fr_average_waiting_5 = [];
% Fr_average_waiting_5 = Fr_changes_test(:,5);
% 
% frM.Fr_average_waiting_1 = Fr_average_waiting_1;
% frM.Fr_average_waiting_2 = Fr_average_waiting_2;
% frM.Fr_average_waiting_3 = Fr_average_waiting_3;
% frM.Fr_average_waiting_4 = Fr_average_waiting_4;
% frM.Fr_average_waiting_5 = Fr_average_waiting_5;


% - - - - - - - -
% - - - Computation of the Depth
% - - - - - - -

for neuron = 1:size(frM,1)
    if isnan(frM.Id_spont(neuron))
        load(fullfile(sort_horridge{frM.Recording(neuron)},'chanMap.mat'));
        clu = frM.Id_horr(neuron);
        idepth = get_cluster_depth(sort_horridge{frM.Recording(neuron)}, clu);
    else
        load(fullfile(sort_before{frM.Recording(neuron)},'chanMap.mat'));  
        clu = frM.Id_spont(neuron);
        idepth = get_cluster_depth(sort_before{frM.Recording(neuron)}, clu);
    end
    if master
        depths_M(neuron) = abs(idepth.depth' - ycoords(masters_depth(frM.Recording(neuron))));
    else
        depths_M(neuron) = abs(idepth.depth' - ycoords(yokes_depth(frM.Recording(neuron))));
    end
end
frM.depth = depths_M';
if master
    save_var = "masters";
else
    save_var = "yokes";
end
% up_regulated
save(sprintf("./memory_cells/height_spike_rate/%i.mat",recording_number),"recording_number","frM","heights")
end
