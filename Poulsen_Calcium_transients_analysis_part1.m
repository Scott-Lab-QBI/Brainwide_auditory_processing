MatFiles=dir('*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
Noise=load(name, 'Noise');
Noise=Noise.Noise;

Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.

GoodNoise=Noise(Fitness,:);

Rs=load(name, 'ROIs');
Rs=Rs.ROIs;Rs=Rs(:,Fitness);
cor_name=strrep(name,'analysis_matlab','correlation');
cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
dims=size(cor_im);
ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
Centroids=zeros(size(Rs,2),2);
for roi_nb=1:size(ROI,3)
    progressbar([],roi_nb/size(ROI,3));
    temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');    
    temp=temp.Centroid;
    Centroids(roi_nb,1:2)=temp;
end
MatFiles(1).ROIs=Centroids;

MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
    progressbar(i/length(MatFiles),[]);
    name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
    C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
    C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1
    GC=C(F,:);
    Noise=vertcat(Noise,N);
    GN=N(F,:);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    GoodNoise=vertcat(GoodNoise,GN);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;Rs=Rs(:,F);
    cor_name=strrep(name,'analysis_matlab','correlation');
    cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
    dims=size(cor_im);
    ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
    Centroids=zeros(size(Rs,2),2);
    for roi_nb=1:size(ROI,3)
        progressbar([],roi_nb/size(ROI,3));
        temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
        temp=temp.Centroid;
        Centroids(roi_nb,1:2)=temp;
    end
    MatFiles(i).ROIs=Centroids;
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore
clearvars Rs ROI roi_nb Centroids cor_im cor_name temp dims Calcium Noise Fitness
%%
ZS=zscore(GoodCalcium,1,2); %%%to normalize the data

save('Reb_AudProf_Raw.mat','GoodCalcium', 'GoodNoise', 'MatFiles','-v7.3'); clearvars GoodCalcium GoodNoise GoodCalNoise Fitness Calcium Noise

spike=[0,1.69644104899772,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike);
RebRegressorAudProf=zeros(12,size(ZS,2));counter=1;
%Media1 1200 timepoints
framerate=4;
White_noise=[62];
for i=1:1:(20*framerate)
    idx=(White_noise*framerate)+i;
    RebRegressorAudProf(1,idx)=i;
end
Full_vol=[110 120 130 140 150	160	170	180	190	200];counter=2;
for i=Full_vol
    idx=i*framerate;idx=idx;
    RebRegressorAudProf(counter,idx-1:idx-1+length(spike)-1)=spike';
    counter=counter+1;
end
Vol_sweep=[230	240	250	260	270	280	290]+10;
for i=Vol_sweep
    idx=i*framerate;idx=idx;
    RebRegressorAudProf(12,idx-1:idx-1+length(spike)-1)=spike';
end
%Media2 1080 timepoints
White_noise=[22];
for i=1:1:(20*framerate)
    idx=(White_noise*framerate)+i+1200;
    RebRegressorAudProf(1,idx)=i;
end
Full_vol=[160 150 140 130 120 110 100 90 80 70];counter=2;
for i=Full_vol
    idx=i*framerate;idx=idx+1200;
    RebRegressorAudProf(counter,idx-1:idx-1+length(spike)-1)=spike';
    counter=counter+1;
end
Vol_sweep=[200 240 190 230 210 220 250];
for i=Vol_sweep
    idx=i*framerate;idx=idx+1200;
    RebRegressorAudProf(12,idx-1:idx-1+length(spike)-1)=spike';
end

%Media3 1560 timepoints
White_noise=[222];
for i=1:1:(20*framerate)
    idx=(White_noise*framerate)+i+1200+1080;
    RebRegressorAudProf(1,idx)=i;
end
Full_vol=[90 30 10 100 60 40 20 70 80 50];counter=2;
for i=Full_vol
    idx=i*framerate;idx=idx+1200+1080;
    RebRegressorAudProf(counter,idx-1:idx-1+length(spike)-1)=spike';
    counter=counter+1;
end
Vol_sweep=[190 fsave180	170	160	150	140	130];
for i=Vol_sweep
    idx=i*framerate;idx=idx+1200+1080;
    RebRegressorAudProf(12,idx-1:idx-1+length(spike)-1)=spike';
end
Frequency_sweep_up=[330];
for i=1:1:(30*framerate)
    idx=(Frequency_sweep_up*framerate)+i+1200+1080;
    RebRegressorAudProf(13,idx)=i;
end
Frequency_sweep_down=[280];
for i=1:1:(30*framerate)
    idx=(Frequency_sweep_down*framerate)+i+1200+1080;
    RebRegressorAudProf(14,idx)=(30*framerate)-i;
end
clearvars counter i GCaMP6 spike framerate
RebRegressorAudProf([1 13 14],:)=RebRegressorAudProf([1 13 14],:)./max(RebRegressorAudProf([1 13 14],:),[],2);

parfor i=1:size(ZS,1)
    mdl=fitlm(RebRegressorAudProf',ZS(i,:));    
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

%% Kmeans

idx_rsq01=find([model_basic.rsquared]>0.1);ZS_rsq=ZS(idx_rsq01,:);
figure;plot(mean(ZS(idx_rsq01,100:3750),1));
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS(idx_rsq01,:),50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,RebRegressorAudProf,idxKmeans_ZS,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

GoodClust_goodmembers=[];Threshold=0.5;
for i=1:length(GoodBetas_ZS)    
    ZS_temp=ZS_rsq(find(idxKmeans_ZS==GoodBetas_ZS(i)),:);
    corr_temp=zeros(1,length(ZS_temp));
    parfor jj=1:size(ZS_temp,1)
        temp_corr=corrcoef(Cmap_ZS(GoodBetas_ZS(i),:), ZS_temp(jj,:));
        corr_temp(jj)=temp_corr(1,2);
    end    
    GoodClust_goodmembers(i).ZS=ZS_temp(find(corr_temp>Threshold),:);
    idx_temp=find(idxKmeans_ZS==GoodBetas_ZS(i));
    GoodClust_goodmembers(i).idx=idx_temp(find(corr_temp>Threshold));
    GoodClust_goodmembers(i).mean=mean(GoodClust_goodmembers(i).ZS,1);
    GoodClust_goodmembers(i).STD=std(GoodClust_goodmembers(i).ZS,1,1);       
end

Numbers=[0 [MatFiles.GoodNumber]];
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one
    name=strcat(MatFiles(i).name);    
    [name2,~]=regexp(name,'AuditoryProfiling_(\d+)_','tokens','match'); %%%to get the number of the fish
    [name3,~]=regexp(name,'fish(\d+)_','tokens','match'); %%%to get the number of the fish    
    Fish=strcat(name2{1}{1},name3{1}{1});Fish=str2double(Fish); %%%to get the number of the fish
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter name2 name3 rows %%%to get rid of vairables we will not use anymore
Fish_list=unique(idx_Fish);

%trace,raster and fish contribution figure
counter=1;
for i=1:length(GoodBetas_ZS)    
    idx_temp=GoodClust_goodmembers(i).idx;
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 1800, 900]);
    subplot(1,3,1);
    plot(mean(ZS_rsq(idx_temp,:),1));ylim([-2 5]);    
    subplot(1,3,2);
    imagesc(ZS_rsq(idx_temp,:));
    subplot(1,3,3);
    histogram(categorical(idx_Fish(idx_temp)));
    %print(Fighandle,strcat('\\D\Rebecca\Auditory Profiling\Results\','Profile_rsq_Clust_',num2str(i)),'-dsvg','-r0');
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

for i=1:length(GoodBetas_ZS)    
    idx_temp=GoodClust_goodmembers(i).idx;
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 1200, 1200]);    
    plot(mean(ZS_rsq(idx_temp,:),1));ylim([-2 5]);  
    %print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Rebecca\Auditory Profiling\Results\','Profile_rsq_Clust_mean_',num2str(i)),'-dpdf','-bestfit','-r0');
    close;
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(model_basic)%%% to make a variable the size of ModelResults
    coef=[model_basic(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            %if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            %end
        end
    end
end

idxempty=cellfun('isempty',coefficients); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients=cell2mat(coefficients); %%%to make a matrix of the coefficients array
mean_coef=mean(coefficients,1);
std_coef=std(coefficients,1,1);

coefficients_rsq=coefficients(idx_rsq01,:);
Freq_clust=zeros(10,size(ZS,2));
for freq_coef=1:10
    Freq_clust(freq_coef,:)=mean(ZS_rsq(coefficients_rsq(:,freq_coef+1)>(mean_coef(freq_coef)+2*std_coef(freq_coef)),:),1);
end

%image of each cluster with a trace, raster and fish contributions
for freq_coef=1:10    
    Fighandle=figure;
    set(Fighandle, 'Position', [100, 100, 1200, 1200]);    
    idx_temp=coefficients_rsq(:,freq_coef+1)>(mean_coef(freq_coef)+2*std_coef(freq_coef));
    subplot(1,3,1);
    plot(mean(ZS_rsq(idx_temp,:),1));ylim([-2 5]);    
    subplot(1,3,2);
    imagesc(ZS_rsq(idx_temp,:));
    subplot(1,3,3);
    histogram(categorical(idx_Fish(idx_temp)));    
    print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Rebecca\Auditory Profiling\Results\','Profile_rsq_FreqSelect_withHist_',num2str(freq_coef)),'-dpng','-r0');
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot


for fish=Fish_list'
    fish_name=num2str(fish);
    IndexC=strfind({MatFiles.name},strcat(fish_name(1:length(fish_name)-1),'_fish',fish_name(end)));
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    Centroids=[];
    for MatfileNb=1:length(MatFiles_fish)
        Plane=MatFiles(MatFiles_fish(MatfileNb)).name;Plane=regexp(Plane,'Slice(\d+)_','tokens');Plane=Plane{1}{1};
        Centroids_temp=MatFiles(MatFiles_fish(MatfileNb)).ROIs;
        Centroids_temp(:,3)=str2double(Plane);
        Centroids=vertcat(Centroids, Centroids_temp);
    end    
    csvwrite(strcat('ROIs_',fish_name,'.csv'),Centroids);
end

CSV_Files=dir('_3Warped*.csv');
ROIs=struct();truth=[];
for i=1:length(CSV_Files)
    temp=csvread(CSV_Files(i).name,1);
    Fishname=regexp(CSV_Files(i).name,'_3Warped(\d+)_fish(\d+).csv','tokens');Fishname=strcat('2018',Fishname{1}{1},Fishname{1}{2});
    ROIs(i).name=str2num(Fishname);
    ROIs(i).coord=temp(:,1:3);
    truth(i)=size(temp,1)==sum(idx_Fish==str2num(Fishname))
end
clearvars i temp Fishname

%% Sorting proper data
temp_nb=0;truth=[];
MatFiles_names={MatFiles.name};
Sorted_ROIs=nan(length(idx_Fish),3);
for fish_nb=1:length(Fish_list)
    temp=num2str(ROIs(fish_nb).name);
    IndexC=strfind({MatFiles.name}, strcat(temp(1:length(temp)-1),'_fish',temp(end)));
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    ROI=ROIs(fish_nb).coord;
    Counter_ROI_coord=1;
    for file_nb=1:length(MatFiles_fish)
        if file_nb>1
            Counter_ROI_coord=Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)));
        end
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];
        end

        Sorted_ROIs(numbersForROIs(1):numbersForROIs(2),:)=ROI(Counter_ROI_coord:Counter_ROI_coord+length(Sorted_ROIs(numbersForROIs(1):numbersForROIs(2)))-1,:);
    end
end
sum(~isnan(Sorted_ROIs(:,1)))==length(idx_Fish)
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

ROI_fish=round(Sorted_ROIs);
ROI_WT(:,1)=round(ROI_fish(:,2));
ROI_WT(:,2)=round(ROI_fish(:,1));
ROI_WT(:,3)=round(ROI_fish(:,3)/1.5);

load('C:\Data\Zbrain_Masks.mat')
Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Stratum'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_WT,Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);    
end
clearvars i j k l m temp Mask Hindbrain_Mask IsInBrainRegion IsInEyes_temp IndexC

rsq_WT=[model_basic.rsquared];
ROI_WT_rsq=ROI_WT(idx_rsq01,:);
rsq_WT_rsq=rsq_WT(idx_rsq01);
mean_coef_rsq=mean(coefficients(idx_rsq01,:),1);
std_coef_rsq=std(coefficients(idx_rsq01,:),1,1);
for freq_coef = 1:10       
    temp_WT=coefficients_rsq(:,freq_coef+1);
    idx_WT=find(temp_WT>(mean_coef_rsq(freq_coef)+2*std_coef_rsq(freq_coef)));
    temp_WT=temp_WT/max(temp_WT);    
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Rebecca\Auditory Profiling\Results\ROI_coef_thresh_rsq_',num2str(freq_coef),'.csv'),[ROI_WT_rsq(idx_WT,:) temp_WT(idx_WT) rsq_WT_rsq(idx_WT)']);
    IsInBrain_WT=find(ismember(ROI_WT_rsq(idx_WT,:),Zbrain_AllMask,'rows'));
    csvwrite(strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Rebecca\Auditory Profiling\Results\ROI_coef_thresh_rsq_',num2str(freq_coef),'_inbrain.csv'),[ROI_WT_rsq(idx_WT(IsInBrain_WT),:) temp_WT(idx_WT(IsInBrain_WT)) rsq_WT_rsq(idx_WT(IsInBrain_WT))']);  
end

%% Plot scatter per fish

figure(1);
for i=1:length(Fish_list)
    scatter(ROI_WT_rsq(idx_Fish_rsq==Fish_list(i),1),ROI_WT_rsq(idx_Fish_rsq==Fish_list(i),2),'.');%hold on;
    axis([0 1400 0 800]);
    pause
end

%% Rebecca chose the following clusters
GoodBetas_select=([1 2 5 6 8 11]);
ClusterMean=zeros(length(GoodBetas_select),size(ZS_rsq,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 800]);
for i=1:length(GoodBetas_select)    
    ClusterMean(i,:)=GoodClust_goodmembers(GoodBetas_select(i)).mean;
end
imagesc(ClusterMean,[-1 4]);colormap hot;