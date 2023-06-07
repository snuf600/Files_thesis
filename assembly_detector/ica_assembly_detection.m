%% FastICA source code: https://github.com/aludnam/MATLAB/tree/master/FastICA_25 
%Largely based on :
% Lopes-dos-Santos, V., Ribeiro, S., &#38; Tort, A. B. L. (2013). Detecting cell assemblies in large neuronal populations. Journal of Neuroscience Methods, 220 (2), 149–166. https://doi.org/10.1016/J.JNEUMETH.2013.04.010

function [predicted_nbr_assemblies, predicted_nbr_neurons,assemblies,activity] = ica_assembly_detection(M,plotter)
%given a matrix A containing the zscores of N neurons throughout M time bins,
% this function will detect neuronal assemblies, and the neurons that are
%part of the different assemblies.

%Pearson correlation coefficient from this matrix A
A=M;
[R, ~,~,~] = corrcoef(A);




%Eigenvalues and eigenvectors
E = eig(R);
[V,D] = eig(R);

q = size(A,1)/size(A,2);

%compute the distribution limits
lambda_max = (1 + sqrt(1/q))^2;
lambda_min = (1 - sqrt(1/q))^2;

%eigenvalues above the max limit predict the number of assemblies, while
%eigenalues above+below predict the number of neurons participating in
%these assemblies
predicted_nbr_assemblies = size(find(E > lambda_max),1)
predicted_nbr_neurons_evs = size(find(E > lambda_max | E < lambda_min),1);


%eigenvectors length computing

neuron_vectors = V(:,size(V,1) + 1-predicted_nbr_assemblies:end );
evs = (D((size(V,1) + 1-predicted_nbr_assemblies:end),(size(V,1) + 1-predicted_nbr_assemblies:end)));


reduced_data = neuron_vectors*neuron_vectors'*A';
[~, M, ~] = fastica(reduced_data,'pcaE',neuron_vectors,'pcaD',evs,'verbose','off');%A contains the importance of each neuron in the independent component (kind of new neuron vector)
M_save = M;


%Important components have the same sign since data is all positive and
%otherwise they would cancel each other out
for i = 1:predicted_nbr_assemblies
    [~, idxs] = maxk(abs(M(:,i)),5);%values of 5 max's to check if these are positive or negative
    if sum(M(idxs,i))<0
        M(M(:,i)>0,i)=0;
    else
        M(M(:,i)<0,i)=0;
    end
end

%Plot points if you want
plot(vecnorm(M,2,2),zeros(1,size(M,1)),"k.")

%Find how many neurons are belonging to assemblies by clustering the
%norm of the different neurons (norm of importance in all components
p_column = reshape(vecnorm(M,2,2),[],1);
[idx,~,~,~] = kmeans(p_column,2, 'Replicates', 10);%can be improved by taking into account that negative clusters will be together

mean_2 = mean(p_column(idx == 2));
mean_1 = mean(p_column(idx == 1));

if abs(mean_2)>abs(mean_1)
    thr = min(p_column(idx==2));
else
    thr = min(p_column(idx==1)); %= max(p_column(idx == 1))+0.01;
end


%Find the most important neurons
neurons_idxs = find(p_column >= thr)';
predicted_nbr_neurons =length(neurons_idxs);%= max(size(find(E > lambda_max | E < lambda_min),1),length(neurons_idxs))%Either highest ica 'evs' or deviations from random distribution of pca evs
% [~,neurons_idxs] = maxk(p_column,predicted_nbr_neurons);
% neurons_idxs = sort(neurons_idxs);



%Scale the biggest component to 1 to overcome problem of multiple neurons
%in 1 assembly, we know for these neurons that they belong to at least one
%assembly so scale to the one that's most likely
M_scaled = M(neurons_idxs,:);%./max(abs(M(neurons_idxs,:)),[],2);



%K-means clustering for thresholding in every column
new_idx = zeros(predicted_nbr_neurons,predicted_nbr_assemblies);
for i=1:predicted_nbr_assemblies
    p_column = M_scaled(:,i);
    
    [idx,~,~,~] = kmeans(p_column,2, 'Replicates', 1000);%can be improved by taking into account that negative clusters will be together
    
    
    mean_2 = mean(p_column(idx == 2));
    mean_1 = mean(p_column(idx == 1));
    
    
    if abs(mean_2)>abs(mean_1)
        new_idx(idx==2,i) = 1;
        new_idx(idx==1,i) = 0;
    else
        new_idx(idx==1,i) = 1;
        new_idx(idx==2,i) = 0;
    end
end
%create the assemblies
assembly_vector = new_idx;%reshape(new_idx,[],predicted_nbr_assemblies);
assemblies = cell(predicted_nbr_assemblies,1);
for j = 1:predicted_nbr_assemblies
    assemblies{j} = neurons_idxs(assembly_vector(:,j)==1);
end

if plotter
    figure
    subplot(2,1,1)
    imagesc(A')
    caxis([0 5])
    subplot(2,1,2)
    icacomp = A*M_save;
    plot(icacomp(:,:))
    xlim([0 size(A,1)])
end
activity = icacomp;
end