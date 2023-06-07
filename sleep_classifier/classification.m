plotting = 1;
tsneing=1;
if ~exist("D","var")

    recording_path = "W:\Jeremy\EEG\recordings\JC230309";

    [~,paths]=dos([strcat("dir /s /b ", fullfile(recording_path,'*.oebin'))]);%https://nl.mathworks.com/matlabcentral/answers/249504-i-want-to-get-search-for-all-of-the-mat-files-in-an-exact-directory-and-all-of-the-subdirectories
    paths=textscan(paths,'%s','delimiter','\n'); paths=paths{:};
    nb_files = length(paths);
    
    eeg = [];
    emg = [];

    for i=1:nb_files
        i
        labels = ["NREM","REM","Wake","Unknown"];
        path = paths{i};
        % Load the files
    
        D=load_open_ephys_binary(path,'continuous',1,'mmap');
        recording_fs = round(1/(D.Timestamps(2)-D.Timestamps(1)))%Hz
        fs = 1000;
        downsampling = recording_fs/fs;
        eeg = [eeg double(downsample(D.Data.Data.mapped(1,:),downsampling))];
        emg = [emg double(downsample(D.Data.Data.mapped(2,:),downsampling)-downsample(D.Data.Data.mapped(3,:),downsampling))];

    end
end
epoch_size = 8;%in seconds
epoch_count = 10800;

freq = 1:floor(fs/2);

eeg_reshaped = reshape(eeg(1:floor(length(eeg)/fs/epoch_size)*fs*epoch_size),epoch_size*fs,[]);
emg_reshaped = reshape(emg(1:floor(length(emg)/fs/epoch_size)*fs*epoch_size),epoch_size*fs,[]);


%Calculate power
eeg_pow = getPower(eeg_reshaped,fs,epoch_size);
emg_pow = getPower(emg_reshaped,fs,epoch_size);

%detect abnormal data
eeg_sum = log10(sum(eeg_pow,1));
eeg_sum = eeg_sum-median(eeg_sum);
emg_sum = log10(sum(emg_pow,1));
emg_sum = emg_sum-median(emg_sum);
good_epochs = eeg_sum<1 & emg_sum<2;

eeg_pow = normalize(log10(eeg_pow(:,good_epochs)),1);
emg_pow = normalize(log10(emg_pow(:,good_epochs)),1);

concat = [eeg_pow;emg_pow]';

[pcs,scores,latents] = pca(concat);

if(plotting)
    thin = 1;
    figure
    plot(scores(1:thin:end,1),scores(1:thin:end,2),"k.")
    title("first 2 principal components")
    xlabel("PC1")
    ylabel("PC2")
end

if (plotting)
    figure
    plot3(scores((1:thin:end),1),scores((1:thin:end),2),scores((1:thin:end),3),"b.");
    hold on
    title("3 PCs with ground truth")
    xlabel("PC1")
    ylabel("PC2")
    zlabel("PC3")
end
if(tsneing)% && ~exist("Y","var"))
    Y =  tsne(concat,"NumDimensions",3,'Perplexity',200);%scores(:,1:5),'NumDimensions',3);%so far scores(:,1:5), 'Perplexity',200
    %plot(Y(:,1),Y(:,2),"x")
end
if plotting
    thin = 1;
    figure
    plot3(Y((1:thin:end),1),Y((1:thin:end),2),Y((1:thin:end),3),"b.");
    title("t-SNE with ground truth")
    xlabel("t-SNE1")
    ylabel("t-SNE2")
    zlabel("t-SNE3")
end

for i = 1:0.1:10
    temp_labels = dbscan(Y,i,10);
    if max(temp_labels) == 2
        break;
    end
end
%REM sleep seperated by t-SNE
REM_label = mode(temp_labels(temp_labels ~= mode(temp_labels)));

%Fit gaussian mixture model to data --> PCA will be normally distributed
%around mean
idx = zeros(size(scores,1),1);
idx(temp_labels==REM_label) = 3;
gmfit = fitgmdist(scores(temp_labels~=REM_label,1:3),2,'Replicates',10,'CovarianceType','full', 'SharedCovariance',false); % Fitted GMM
idx(temp_labels~=REM_label) = cluster(gmfit,scores(temp_labels~=REM_label,1:3));

  
function pow = getPower(v,fs,epoch_size)
        %detrendp
%     detrend(v);
%     %hanning
%     v = v.*hann(fs*epoch_size);
    %calculate power
    xdft = fft(v);
    N = size(v,1);
    xdft = xdft(1:N/2+1,:);
    psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx(2:end-1,:) = 2*psdx(2:end-1,:);
    pow = psdx;
    nb_points = fs*epoch_size;
    pow = pwelch(v,100,50,nb_points,fs);
    %drop the first 0.5Hz components
    pow = pow(round(0.5/(fs/2)*(nb_points/2+1)):end,:);
end