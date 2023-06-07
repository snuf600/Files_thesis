clear

recording_database_updated_20200805_JCfreshStart



% i = number of the recording
i = 6;
master = 0;%masters or yokes
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
for i = 1:length(sort_horridge)
    clearvars -except i sort_horridge master sort_before
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
        lim_start=1;
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
    %%
    % Create columm :
    % Recording number
    % Id horridge / Id spont / Id after (same on long recording)
    L = 0.005; %change to increase/decrease time bins
    T = 0.1;%time to inspect before stimulation
    T_wait = 0.0;%time to wait before stimulation so no contamination of stimulation in data
    step = 0.1;
    cont = 0;
    spike_spacing = [];
    long_pause = [];
    
    
    temp1 = events.onsets;%look at all stims to determine waiting time
    temp1 = temp1(tbin(2,1)<=temp1 & temp1<tbin(4,2));
    temp2 = events.offsets;
    temp2 = temp2(tbin(2,1)<=temp2 & temp2<tbin(4,2));
    if ~(min(temp1(2:end)-temp2(1:end-1))>T+T_wait)
        disp("Stimulations not 50ms spaced apart")%check that the stimulation do not interfere
        T = floor((min(temp1(2:end)-temp2(1:end-1))-T_wait)/L)*L;
    end

    for j = 1:length(clusters)
        cont = cont+1;
        if master
            disp(sprintf('master %d, cluster long %d',i,j))
        else
            disp(sprintf('yoke %d, cluster long %d',i,j))
        end
        Recording(cont) = i;
        Id_horr(cont) = clusters(j);
        Id_spont(cont) = clusters(j);
        Id_after(cont) = clusters(j);
        
        % Create the FR column for each phases d epending on Tbin
        for k = 2:4   
            

            stim_ons = events.onsets;
            stim_ons = stim_ons(tbin(k,1)<=stim_ons & stim_ons<tbin(k,2));
            stim_ofs = events.offsets;
            stim_ofs = stim_ofs(tbin(k,1)<=stim_ofs & stim_ofs<tbin(k,2));
            idxs = (stim_ons(2:end)-stim_ofs(1:end-1)> 0.2);
            
            spikes = zeros(length(stim_ofs),round(T/L));
            if length(spike_spacing) <= length(stim_ofs)
                spike_spacing = [spike_spacing min(stim_ons(2:end)-stim_ofs(1:end-1)) reshape(stim_ons(2:end)'-stim_ofs(1:end-1)',1,[])];
                long_pause = spike_spacing > 0.2;
            end
            for stop_stim = 1:length(stim_ons)
                edges = stim_ons(stop_stim)-T_wait-T:L:stim_ons(stop_stim)-T_wait;
                spikes(stop_stim,:) = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
            end

%             while win(2)<=tbin(k,2)
                
                
%                 stj = spk.st(spk.clu == clusters(j) & spk.st>win(1) & spk.st<=win(2));
%                 frmov = [frmov numel(stj)/L];
%                 win = win+L;
%             end
            switch k
%                 case 1
%                     Fr_spont{cont} = frmov;
                case 2
                    Fr_early{cont} = spikes;
                case 3
                    Fr_middle{cont} = spikes;
                case 4
                    Fr_late{cont} = spikes;
%                 case 5
%                     Fr_waiting{cont} = frmov;
            end

            
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
frM.Fr_early = Fr_early';
frM.Fr_middle = Fr_middle';
frM.Fr_late = Fr_late';
frM.Id_after = Id_after';
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
up_down_regulation
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for i = 1:length(sort_horridge)%; 5 is 23/09
    clearvars -except i sort_horridge master sort_before spk ycoords masters_depth yokes_depth
    recording_number = i;


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
%     try
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
        
%         %Load first and last event for Eearly and Late
%         %For it: look at the event file.
%         if i==1 %25/07/22 // Ptf1a
%             lim_start=1;
%             lim_event=324;
%     
%         elseif i==2 %19/08/22 // Ptf1a
%             lim_start=1;
%             lim_event=336;
%         
%         elseif i==3 %19/08/22 // Eng (waiting: 3676-7264)
%             lim_start=1;
%             lim_event=749;
%             waiting_start = 3676;
%             waiting_stop = 7264;
%         elseif i==4 %21/09/22 // Ptf1a (waiting: 3709-7382)
%             lim_start=2;
%             lim_event=1070;
%             waiting_start = 3709;
%             waiting_stop = 7382;
%         elseif i==5 %23/09/22 // Ptf1a 
%             lim_start=1;
%             lim_event=128;
%             waiting_start = 3921;
%             waiting_stop = 7924;
%         elseif i==6 %11/10/22 // Ptf1a 
%             lim_start=1;
%             lim_event=103;
%             waiting_start = 3940;
%             waiting_stop = 7563;
%         elseif i==7 %16/11/22 // Ptf1a 
%             lim_start=3;
%             lim_event=419;
%             waiting_start = 3025;
%             waiting_stop = 6902;
%         elseif i==8 %24/11/22 // ENG 
%             lim_start=1;
%             lim_event=164;
%             waiting_start = 3075;
%             waiting_stop = 7212;
%         elseif i==9 %28/11/22 // ENG
%             lim_start=1;
%             lim_event=225;
%             waiting_start = 3207;
%             waiting_stop = 7213;
            
            
            
            
            % Initialization parameter in case of problem
%         else
            lim_start=1;
            
            while events(evtid).onsets(lim_start+1)>events(evtid).onsets(lim_start)+1 || events.onsets(lim_start)<900
                lim_start = lim_start+1;
                if lim_start == 100
                    break
                end
            end
            assert(events.onsets(lim_start)>=900)
            lim_event=sum(events(evtid).onsets(lim_start) <= events(evtid).onsets & events(evtid).onsets <= events(evtid).onsets(lim_start)+600)+lim_start-1;
%         end
        
        
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
%         disp(tbin);
        
        %%
        % Create columm :
        % Recording number
        % Id horridge / Id spont / Id after (same on long recording)
        L = 0.005; %change to increase/decrease time bins
        T = 0.1;%time to inspect before stimulation
        T_wait = 0.0;%time to wait before stimulation so no contamination of stimulation in data
        step = 0.1;
        cont = 0;
        spike_spacing = [];
        long_pause = [];
        
        
        temp1 = events.onsets;%look at all stims to determine waiting time
        temp1 = temp1(tbin(2,1)<=temp1 & temp1<tbin(4,2));
        temp2 = events.offsets;
        temp2 = temp2(tbin(2,1)<=temp2 & temp2<tbin(4,2));
        if ~(min(temp1(2:end)-temp2(1:end-1))>T+T_wait)
            disp("Stimulations not 50ms spaced apart")%check that the stimulation do not interfere
            T = floor((min(temp1(2:end)-temp2(1:end-1))-T_wait)/L)*L;
        end
    
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
            for k = 2:4   
                
    
                stim_ons = events.onsets;
                stim_ons = stim_ons(tbin(k,1)<=stim_ons & stim_ons<tbin(k,2));
                stim_ofs = events.offsets;
                stim_ofs = stim_ofs(tbin(k,1)<=stim_ofs & stim_ofs<tbin(k,2));
%                 idxs = (stim_ons(2:end)-stim_ofs(1:end-1)> 0.2);
                
                spikes = zeros(length(stim_ofs),round((2*T)/L));
%                 if length(spike_spacing) <= length(stim_ofs)
%                     spike_spacing = [spike_spacing min(stim_ons(2:end)-stim_ofs(1:end-1)) reshape(stim_ons(2:end)'-stim_ofs(1:end-1)',1,[])];
%                     long_pause = spike_spacing > 0.2;
%                 end
                for stop_stim = 1:length(stim_ons)
                    edges = stim_ons(stop_stim)-T_wait-T:L:stim_ons(stop_stim)-T_wait+T;
                    spikes(stop_stim,:) = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
                end
    
  
                switch k
   
                    case 2
                        Fr_early{cont} = spikes;
                        edges = tbin(k,1):L:tbin(k+2,2);
                        frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
                        Fr_horridge{cont} = frmov;
                    case 3
                        Fr_middle{cont} = spikes;
                    case 4
                        Fr_late{cont} = spikes;
   
                end
    
                
            end
%             k = 1;
%             edges = tbin(k,1):L:tbin(k,2);
%             frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
%             Fr_spont{cont} = frmov;
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
            
%             Fr_average_spont(cont)=mean(Fr_spont{cont});
            
    %         Fr_average_waiting(cont) = mean(Fr_waiting{cont});
    
            Tbin_tot{cont}=tbin;
        end
    
    
    
    %%
    % Create the frM table
    
    
    frM = table();
    frM.Recording = Recording';
    frM.Tbin = Tbin_tot';
    frM.Id_spont = Id_spont';
%     frM.Fr_spont = Fr_spont';
%     frM.Fr_average_spont = Fr_average_spont';
    frM.Id_horr = Id_horr';
    frM.Fr_early = Fr_early';
    frM.Fr_middle = Fr_middle';
    frM.Fr_late = Fr_late';
    frM.Id_after = Id_after';
    frM.Fr_Horridge = Fr_horridge';
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
    
    
%     - - - - - - - -
%     - - - Computation of the Depth
%     - - - - - - -
%     visualise_spikes_after_stim
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
    load(sprintf('D:/From c disk/Documenten/Arnaud/School/5e jaar/NERF/files train/workspaces/beforevafter/upregulatedxd/master/%i/data.mat',i))
    frM.depth = depths_M';
    save(sprintf('D:/From c disk/Documenten/Arnaud/School/5e jaar/NERF/files train/workspaces/beforevafter/upregulatedxd/master/%i/data.mat',i),"L","events","frM","p2","p_value","relevant_neurons","stim_number")


%     up_regulated
%         
%     downregulated_after_stim
%     catch
% %         continue
%     end
end


% 
%% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

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
upup = 0;
%%
for i = 1:length(sort_horridge)%; 5 is 23/09
    clearvars -except i sort_horridge master sort_before upup
    try
        i
        recording_number = i;
        if master
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
        else
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
        end
        downregulated_after_stim
    catch
        continue
    end
end
%%
recording_database_updated_20200805_JCfreshStart

% i = number of the recording
i = 5;
master = 0;%masters or yokes
if master
    sort_horridge = sort_masters_horridge;
    sort_before = sort_masters_before;
else
    sort_horridge = sort_yokes_horridge;
    sort_before = sort_yokes_before;
end
upup = 0;

%%
for i = 1:length(sort_horridge)%; 5 is 23/09
    clearvars -except i sort_horridge master sort_before upup
    try
        recording_number = i
        if master
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
        else
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
        end
        downregulated_after_stim
    catch
        continue
    end
end


%%
recording_database_updated_20200805_JCfreshStart
upup = 1;
% i = number of the recording
i = 5;
master = 0;%masters or yokes
if master
    sort_horridge = sort_masters_horridge;
    sort_before = sort_masters_before;
else
    sort_horridge = sort_yokes_horridge;
    sort_before = sort_yokes_before;
end

%%
for i = 1:length(sort_horridge)%; 5 is 23/09
    clearvars -except i sort_horridge master sort_before upup
    try
        recording_number = i
        if master
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
        else
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
        end
        downregulated_after_stim
    catch
        continue
    end
end



%%
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
upup = 1
%%
for i = 1:length(sort_horridge)%; 5 is 23/09
    clearvars -except i sort_horridge master sort_before upup
    try
        recording_number = i
        if master
            load(sprintf("workspaces/beforevafter/upregulatedxd/master/%i/data.mat",i))
        else
            load(sprintf("workspaces/beforevafter/upregulatedxd/yoke/%i/data.mat",i))
        end
        downregulated_after_stim
    catch
        continue
    end
end




