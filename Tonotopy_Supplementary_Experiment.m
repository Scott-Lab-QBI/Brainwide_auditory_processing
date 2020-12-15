MatFiles=dir('*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
%MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
Rs=load(name, 'ROIs');
Rs=Rs.ROIs;Rs=Rs(:,Fitness);
cor_name=strrep(name,'analysis_matlab','correlation');
cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
dims=size(cor_im);
ROI_temp=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
Centroids=zeros(size(Rs,2),2);
for i=1:size(ROI_temp,3)
    temp=squeeze(ROI_temp(:,:,i));   
    s=regionprops(temp>0,temp,'WeightedCentroid');
    Centroids(i,:)=s.WeightedCentroid;    
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
    ROI_temp=reshape(full(Rs),dims(1),dims(2),size(Rs,2));
    Centroids=zeros(size(Rs,2),2);
    for ij=1:size(ROI_temp,3)
        progressbar([],ij/size(ROI_temp,3));
        temp=squeeze(ROI_temp(:,:,ij));
        s=regionprops(temp>0,temp,'WeightedCentroid');
        Centroids(ij,:)=s.WeightedCentroid;
    end
    MatFiles(i).ROIs=Centroids;
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore
clearvars Rs ROI roi_nb Centroids cor_im cor_name temp dims Calcium Noise Fitness

%%

%GoodCalNoise=GoodCalcium+GoodNoise;

% ZS=zscore(GoodCalcium,1,2); %%%to normalize the data
% ZS=detrend(ZS')';%%% to Remove a linear trend from ZS 
 
ZS=zscore(GoodCalcium,1,2); %%%to normalize the data

save('Reb_AudProf_Tonotopy_Raw.mat','GoodCalcium', 'GoodNoise', 'MatFiles','-v7.3'); clearvars GoodCalcium GoodNoise GoodCalNoise Fitness Calcium Noise

spike=[0,1.69644104899772,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike);

freqs=sort(unique(Stimuli(:,2)));
RebRegressorAudProf=zeros(length(freqs)+1,size(ZS,2));

for ij=1:size(Stimuli,1)
    idx=Stimuli(ij,1);
    freq_nb=find(freqs==Stimuli(ij,2));
    RebRegressorAudProf(freq_nb,idx-1:idx-1+length(spike)-1)=spike';
    RebRegressorAudProf(length(freqs)+1,idx-1:idx-1+length(spike)-1)=1/Stimuli(ij,3);
end
clearvars counter i GCaMP6 spike framerate

parfor i=1:size(ZS,1)
    mdl=fitlm(RebRegressorAudProf',ZS(i,:));    
    model_basic(i).coef=mdl.Coefficients;        
    model_basic(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i

parfor i=1:size(ZS,1)
    mdl=fitlm(RebRegressorAudProf',ZS(i,:),'interactions');    
    model_interact(i).coef=mdl.Coefficients;        
    model_interact(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars slow_fast slow_fast_fwd i
RebRegressorAudProf(RebRegressorAudProf<0)=0;

rsq_basic=[model_basic.rsquare];
rsq_interact=[model_interact.rsquared];
figure;histogram(rsq_basic);
idx_rsq_interact=rsq_interact>0.2;
idx_rsq_basic=rsq_basic>0.2;

figure;plot(mean(ZS(idx_rsq_interact,100:end),1));hold on;
plot(mean(ZS(idx_rsq_basic,100:end),1));
options = statset('UseParallel',1); [idxKmeans_ZS_interact Cmap_ZS_interact]=kmeans(ZS(idx_rsq_interact,:),50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS_interact,GoodBetas_ZS_interact]=Test_Regress(Cmap_ZS_interact,RebRegressorAudProf,idxKmeans_ZS_interact,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3


options = statset('UseParallel',1); [idxKmeans_ZS_basic Cmap_ZS_basic]=kmeans(ZS(idx_rsq_basic,:),50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS_basic,GoodBetas_ZS_basic]=Test_Regress(Cmap_ZS_basic,RebRegressorAudProf,idxKmeans_ZS_basic,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

%% Linreg selection
coefficients_basic={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(model_basic)%%% to make a variable the size of ModelResults
    coef=[model_basic(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            %if coef.pValue(coef_idx)<0_basic%%%to select the coef that are bellow the p value we want, in this case 0_basic
                coefficients_basic{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            %end
        end
    end
end
idxempty=cellfun('isempty',coefficients_basic); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients_basic(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients_basic=cell2mat(coefficients_basic); %%%to make a matrix of the coefficients array


ZS_rsq=ZS(idx_rsq_basic,:);
mean_coef = mean(coefficients_basic(:,1:10),1);
std_coef = std(coefficients_basic(:,1:10),1,1);

coefficients_rsq_brainReg=coefficients_basic(idx_rsq_basic,:);
sum_Rois2 = zeros(1,4);
std_thresh3=2.5;
std_thresh_exclude=2;
Freq_data25=struct();
for freq_coef=1:10        
    idx_freq=[1:1:10];
    idx_freq(freq_coef)=[];
    %idx_temp=coefficients_rsq_brainReg(:,freq_coef)>(mean_coef(freq_coef)+std_thresh*std_coef(freq_coef));
    idx_temp=(coefficients_rsq_brainReg(:,freq_coef)>(mean_coef(freq_coef)+std_thresh3*std_coef(freq_coef))) ... 
              & all(coefficients_rsq_brainReg(:,idx_freq)<(mean_coef(idx_freq)+std_thresh_exclude*std_coef(idx_freq)),2);    
    Freq_data25(freq_coef).idx=idx_temp;
    Freq_data25(freq_coef).meanZS=mean(ZS_rsq(idx_temp,1:end),1);
    Freq_data25(freq_coef).Raster=ZS_rsq(idx_temp,1:end);
    %Freq_data25(freq_coef).ROIs=ROI_WT_rsq(idx_temp,:);
    %Freq_data(freq_coef).hist=histogram(categorical(idx_Fish(idx_temp)));    
    %print(Fighandle,strcat('\\shares.rdm.uq.edu.au\scottlab-q0291\Rebecca\Auditory Profiling\Results\','Profile_rsq_FreqSelect_withHist_',num2str(freq_coef)),'-dsvg','-r0');
    sum_Rois2(1) = sum_Rois2(1) + sum(idx_temp);
end

%WB trace for all included freq
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10 
    axes(ha(freq_coef));
    plot(Freq_data25(freq_coef).meanZS/max(Freq_data25(freq_coef).meanZS),'Color',color_temp(freq_coef,:),'LineWidth',2);hold on;
    title(strcat('freq :',num2str(freqs(freq_coef)),' - ',num2str(sum(Freq_data25(freq_coef).idx)))); ylim([-0.5 1.1]);xlim([0 1260]);
    for ij=1:length(Stimuli)
        idx=Stimuli(ij,1);
        freq_nb=find(freqs==Stimuli(ij,2));
        if freq_nb==freq_coef
            rectangle('Position',[Stimuli(ij,1) -0.5 8 0.2/Stimuli(ij,3)],'FaceColor',color_temp(freq_nb,:),'EdgeColor',color_temp(freq_nb,:));
        end
    end
end
%print(Fighandle,strcat('..\Figures\','mean_zs_per_tone_wholebrain_crit1_100-end.svg'),'-dsvg','-r0');

%% Creates index

Numbers=[0 [MatFiles.GoodNumber]];
idx_Fish=nan(length(ZS),1);
for i=1:length(MatFiles) %%%%to take slices one by one
    name=strcat(MatFiles(i).name);    
    [name2,~]=regexp(name,'Tonotopy_2019(\d+)_fish(\d+)_','tokens'); %%%to get the number of the fish    
    Fish=strcat(name2{1}{1},name2{1}{2});Fish=str2double(Fish); %%%to get the number of the fish
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter name2 name3 rows %%%to get rid of vairables we will not use anymore
Fish_list=unique(idx_Fish);


%% Prep ROIs for warp
for fish=Fish_list'
    fish_name=num2str(fish);
    IndexC=strfind({MatFiles.name},strcat(fish_name(1:length(fish_name)-2),'_fish',fish_name(end-1:end)));
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

%% Load the warped ROIs
% 1202_2 is crap, 1202_7 is not great, 1202_3 is very tilted
CSV_Files=dir('_2Warped*.csv');
ROIs=struct();truth=[];
for i=1:length(CSV_Files)
    temp=csvread(CSV_Files(i).name,1);
    Fishname=regexp(CSV_Files(i).name,'_2Warped(\d+).csv','tokens');Fishname2=Fishname{1}{1};Fishname=strcat('2019',Fishname2(1:end-2),'_fish',Fishname2(end-1:end));
    ROIs(i).name=str2num(Fishname2);
    ROIs(i).coord=temp(:,1:3);
    truth(i)=size(temp,1)==sum(idx_Fish==str2num(Fishname2))
end
clearvars i temp Fishname

%% Sorting proper data
temp_nb=0;truth=[];
MatFiles_names={MatFiles.name};
Sorted_ROIs=nan(length(idx_Fish),3);
for fish_nb=1:length(Fish_list)
    temp=num2str(ROIs(fish_nb).name);
    IndexC=strfind({MatFiles.name}, strcat(temp(1:end-2),'_fish',temp(end-1:end)));
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

%% check warp
figure(1);
for i=1:length(Fish_list)
    scatter(ROI_WT(idx_Fish==Fish_list(i),1),ROI_WT(idx_Fish==Fish_list(i),2),'.');%hold on;
    axis([0 1400 0 800]);
    pause
end

%%
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


%% 
ROI_WT_rsq=ROI_WT(idx_rsq_basic,:);


Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1400, 800]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10 
    axes(ha(freq_coef));
    scatter(ROI_WT_rsq(Freq_data25(freq_coef).idx,1),ROI_WT_rsq(Freq_data25(freq_coef).idx,2),'.');
    axis([400 1400 0 700]);
end
%print(Fighandle,strcat('..\Figures\','mean_zs_per_tone_wholebrain_crit1_100-end.svg'),'-dsvg','-r0');


%%
%WB trace for all included freq

idx_Fish_rsq=idx_Fish(idx_rsq_basic);
Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=5;
yplot=2;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for freq_coef = 1:10 
    axes(ha(freq_coef));
    histogram(categorical(idx_Fish_rsq(Freq_data25(freq_coef).idx)),categorical(Fish_list),'normalization','probability');    
    title(strcat('freq :',num2str(freqs(freq_coef)),' - ',num2str(sum(Freq_data25(freq_coef).idx))));    
end


