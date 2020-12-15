xplot=5;
yplot=2;
select_volume=2;threshold_inclusion=3;threshold_exclusion=3;
for select_volume=1:4
    Fighandle=figure;
    set(Fighandle, 'Position', [10,10, 1200, 1200]);
    ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
    for freq_coef = 1:10
        axes(ha(freq_coef));
        idx_freq=[1:1:10];
        idx_freq(freq_coef)=[];
        idx_temp=Max_responses(:,freq_coef,select_volume)>threshold_inclusion & max(Max_responses(:,idx_freq,1),[],2)<threshold_exclusion;        
        histogram(categorical(idx_Fish_rsq(idx_temp)),categorical(Fish_list),'normalization','probability'); 
        title(num2str(sum(idx_temp)));
    end
end

for select_volume=1:4
    Fighandle=figure;
    set(Fighandle, 'Position', [10,10, 1200, 1200]);
    ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
    for freq_coef = 1:10
        axes(ha(freq_coef));
        idx_freq=[1:1:10];
        idx_freq(freq_coef)=[];
        idx_temp=Max_responses(:,freq_coef,select_volume)>threshold_inclusion & max(Max_responses(:,idx_freq,1),[],2)<threshold_exclusion;
        imagesc(squeeze(mean(Max_responses(idx_temp,:,:),1)));
        title(num2str(sum(idx_temp)));
    end
end


Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    idx_freq=[1:1:10];
    idx_freq(freq_coef)=[];
        for select_volume=3:4
        idx_temp=Max_responses(:,freq_coef,select_volume)>threshold_inclusion & max(Max_responses(:,idx_freq,1),[],2)<threshold_exclusion;
        scatter(ROI_WT_rsq(idx_temp,1),ROI_WT_rsq(idx_temp,2),'.'); hold on;
        title(num2str(sum(idx_temp)));
    end
    axis([400 1400 0 700]);
end

%% 
Max_responses_flat=reshape(Max_responses,size(Max_responses,1),size(Max_responses,2)*size(Max_responses,3));
figure;imagesc(Max_responses_flat);

[W,H,D] = nnmf(Max_responses_flat,20,'replicates',5);
figure;imagesc(H);

index_NMF_sort=zeros(1,10);
for i=1:10
    [B,index]=sortrows(H,i:10:i+30);
    index_NMF_sort(i)=index(end);
end
figure;imagesc(H(index_NMF_sort,:));

Cluster_NMF=zeros(10,size(ZS_rsq,2));
mean_W=mean(W,1);
std_WT=std(W,1,1);
idx_NMF={};
for i=1:10
    idx_NMF{i}=find(W(:,index_NMF_sort(i))>(mean_W(index_NMF_sort(i))+2*std_WT(index_NMF_sort(i))));    
    Cluster_NMF(i,:)=mean(ZS_rsq(idx_NMF{i},:),1);
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=1;
yplot=2;
select_volume=1;threshold_inclusion=2;threshold_exclusion=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
axes(ha(1));imagesc(H(index_NMF_sort,:));
axes(ha(2));imagesc(Cluster_NMF);

%%
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1800, 600]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    color_temp=repmat(W(:,index_NMF_sort(freq_coef)),1,3);
    color_temp=color_temp/max(color_temp(:));
    scatter(ROI_WT_rsq(:,1),ROI_WT_rsq(:,2),10,color_temp,'.'); hold on;    
    axis([400 1400 0 700]);
    set(gca,'Color','k')
end


%%

[coeff,score,latent,tsquared,explained,mu] = pca(Max_responses_flat);
figure;imagesc(coeff);

index_PCA_sort=zeros(1,10);
for i=1:10
    [B,index]=sortrows(coeff',i:10:i+30);
    index_PCA_sort(i)=index(end);
end
figure;imagesc(coeff(:,index_PCA_sort));

Cluster_PCA=zeros(10,size(ZS_rsq,2));
mean_score=mean(score,1);
std_score=std(score,1,1);
idx_PCA={};
for i=1:10
    idx_PCA{i}=find(score(:,index_PCA_sort(i))>(mean_score(index_PCA_sort(i))+2*std_score(index_PCA_sort(i))));    
    Cluster_PCA(i,:)=mean(ZS_rsq(idx_PCA{i},:),1);
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=1;
yplot=2;
select_volume=1;threshold_inclusion=2;threshold_exclusion=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
axes(ha(1));imagesc(coeff(:,index_PCA_sort)');
axes(ha(2));imagesc(Cluster_PCA);


%% 

options = statset('UseParallel',1); 
[idxKmeans_MaxResp Cmap_MaxResp]=kmeans(Max_responses_flat,40,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_MaxResp);
index_kmeans_sort=zeros(1,10);
for i=1:10
    [B,index]=sortrows(Cmap_MaxResp,i:10:i+30);
    index_kmeans_sort(i)=index(end);
end
figure;imagesc(Cmap_MaxResp(index_kmeans_sort,:));

Cluster_kmeans=[];%zeros(10,size(ZS_rsq,2));
idx_kmeans={};
for i=1:10
    idx_kmeans{i}=find(idxKmeans_MaxResp==index_kmeans_sort(i));    
    Cluster_kmeans(i,:)=mean(ZS_rsq(idx_kmeans{i},:),1);
end
figure;imagesc(Cluster_kmeans);

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    imagesc(ZS_rsq(idx_kmeans{freq_coef},:),[-1 4]);colormap hot;    
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    histogram(categorical(idx_Fish_rsq(idx_kmeans{freq_coef})),categorical(Fish_list));
end

%% Looks like individual fish, so let's normalize the max response, per fish

options = statset('UseParallel',1); 
Max_responses_flat2=zscore(Max_responses_flat,1,2);
for fish=Fish_list'
    Max_responses_flat2(idx_Fish_rsq==fish,:)=Max_responses_flat2(idx_Fish_rsq==fish,:)./max(
end


[idxKmeans_MaxResp Cmap_MaxResp]=kmeans(Max_responses_flat2,20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
figure;imagesc(Cmap_MaxResp);
index_kmeans_sort=zeros(1,10);
for i=1:10
    [B,index]=sortrows(Cmap_MaxResp,i:10:i+30);
    index_kmeans_sort(i)=index(end);
end
Cluster_kmeans=[];%zeros(10,size(ZS_rsq,2));
idx_kmeans={};
for i=1:10
    idx_kmeans{i}=find(idxKmeans_MaxResp==index_kmeans_sort(i));    
    Cluster_kmeans(i,:)=mean(ZS_rsq(idx_kmeans{i},:),1);
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    imagesc(ZS_rsq(idx_kmeans{freq_coef},:),[-1 4]);colormap hot;    
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    histogram(categorical(idx_Fish_rsq(idx_kmeans{freq_coef})),categorical(Fish_list));
end



%%

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));    
    scatter(ROI_WT_rsq(idx_temp,1),ROI_WT_rsq(idx_temp,2),'.'); hold on;
    title(num2str(sum(idx_temp)));    
    axis([400 1400 0 700]);
end

%% Try to cluster per fish

options = statset('UseParallel',1); 

idxKmeans_MaxResp_fish=cell(length(Fish_list),3);
Cmap_MaxResp_fish=cell(length(Fish_list),3);
idx_KmeansPool=[];
for fish=1:length(Fish_list)
    idx_fish_temp=find(idx_Fish_rsq==Fish_list(fish));
    [idxKmeans_MaxResp_fish{fish,1} Cmap_MaxResp_fish{fish,1}]=kmeans(Max_responses_flat(idx_fish_temp,:),20,'Options',options,'Replicates',3,'MaxIter',1000,'Display','final');
    index_kmeans_sort=zeros(1,10);
    for i=1:10
        [B,index]=sortrows(Cmap_MaxResp_fish{fish,1},i:10:i+30);
        index_kmeans_sort(i)=index(end);
    end   
    idxKmeans_MaxResp_fish{fish,2}=index_kmeans_sort;
    Cluster_kmeans=zeros(10,size(ZS_rsq,2));    
    for i=1:10
        idx_kmeans{i}=find(idxKmeans_MaxResp_fish{fish,1}==index_kmeans_sort(i));
        Cluster_kmeans(i,:)=mean(ZS_rsq(idx_fish_temp(idx_kmeans{i}),:),1);
        idx_KmeansPool=vertcat(idx_KmeansPool,[idx_fish_temp(idx_kmeans{i}) i*ones(size(idx_fish_temp(idx_kmeans{i})))]);
    end
    Cmap_MaxResp_fish{fish,2}=Cluster_kmeans;    
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    idx_temp=idx_KmeansPool(idx_KmeansPool(:,2)==freq_coef,1);
    %imagesc(ZS_rsq(idx_temp(randperm(length(idx_temp))),:),[-1 4]);colormap hot;
    imagesc(ZS_rsq(idx_temp,:),[-1 4]);colormap hot;
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    plot(mean(ZS_rsq(idx_KmeansPool(idx_KmeansPool(:,2)==freq_coef,1),:),1))
    for ij=1:length(Stimuli)
        idx=Stimuli(ij,1);
        freq_nb=find(freqs==Stimuli(ij,2));
        if freq_nb==freq_coef
            rectangle('Position',[Stimuli(ij,1) -0.5 8 0.2/Stimuli(ij,3)],'FaceColor',color_temp(freq_nb,:),'EdgeColor',color_temp(freq_nb,:));
        end
    end
    xlim([0 1260]);
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    idx_temp=idx_KmeansPool(idx_KmeansPool(:,2)==freq_coef,1);    
    histogram(categorical(idx_Fish_rsq(idx_temp)),categorical(Fish_list));    
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1800, 600]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10
    axes(ha(freq_coef));
    idx_temp=idx_KmeansPool(idx_KmeansPool(:,2)==freq_coef,1);    
    scatter(ROI_WT_rsq(idx_temp,1),ROI_WT_rsq(idx_temp,3),40,color_temp(freq_coef,:),'.'); hold on;    
    axis([400 1400 0 200]);
    set(gca,'Color','k')
end

MaxResp_perFish=zeros(10,40,length(Fish_list));
for fish=1:length(Fish_list)
    idx_fish_temp=find(idx_Fish_rsq==Fish_list(fish));
    Max_temp=Cmap_MaxResp_fish{fish,1};
    MaxResp_perFish(:,:,fish)=Max_temp(idxKmeans_MaxResp_fish{fish,2},:);
end
figure;
imagesc(squeeze(median(MaxResp_perFish,3)));


PrismTemp=zeros(40,100);
for i=1:10
    PrismTemp(:,1+10*(i-1):10+10*(i-1))=squeeze(MaxResp_perFish(i,:,:))';
end

