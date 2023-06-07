% onsets = events.onsets(events.onsets<1600 & events.onsets>900)-900
% t = 0:1/30:18000/30-1/30; 
% 
% figure;hold on
% plot(t,ZMaster(:,1),'b')
% plot((onsets),ones(1,length(onsets)),"r*")
% s1 = 0.333333;
% s2 = 499.967;
% t1 = s2-s1;
% t2 = onsets(end)-onsets(2);
% t = t*t2/t1-(s1*t2/t1-onsets(2));
% 
% figure;hold on
% plot(t,ZMaster(:,1),'b')
% plot((onsets),ones(1,length(onsets)),"r*")


[matched,~,~] = find_drop(events,t,zData);

fig = figure
subplot(2,1,1)
zData = ZMaster(:,1);
plot(t-min(t),zData)
xlim([min(t) max(t)])
xlabel(["Time [s]"])
ylabel(["Height [mm]"])
title("Foot height recording of a mouse")
hold on;
xlim([0 max(t)-min(t)])
subplot(2,1,2)
plot(t-min(t),zData)
xlabel(["Time [s]"])
ylabel(["Height [mm]"])
title("Zoom on transition from early to late learning phase")
hold on;
xlim([0 max(t)-min(t)])
plot(t(matched)-min(t),zData(matched),"rx")
legend("Height data", "Stimulations")