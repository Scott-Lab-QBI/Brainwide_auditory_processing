%% INPUTS
% CHANGE THE NEXT LINE TO CHOOSE RSQ THRSHOLD VALUE 
idx_rsq = find([model_basic.rsquared]>0.1);
ZS_rsq = ZS(idx_rsq,:);
mean_coef = mean(coefficients,1);
std_coef = std(coefficients,1,1);

%% Run analysis for whole brain data CRITERIA 1
% Filter ROIs that respond selectively to a single pure tone (> MEAN
% + 2.5*STD) and do not respond to the other tones (< MEAN 
% + 2*STD)
% select only*auditory responsive ROIs (>0.1 r squared)
coefficients_rsq_brainReg = coefficients(idx_rsq,:);
sum_Rois2 = zeros(1,4);
std_thresh3 = 2.5;
std_thresh_exclude = 2;
Freq_data25 = struct();
for freq_coef=1:10        
    idx_freq = [2:1:11];
    idx_freq(freq_coef)=[];
    %idx_temp=coefficients_rsq(:,freq_coef+1)>(mean_coef(freq_coef)+std_thresh*std_coef(freq_coef));
    idx_temp=(coefficients_rsq_brainReg(:,freq_coef+1)>(mean_coef(freq_coef+1)+std_thresh3*std_coef(freq_coef+1))) ... 
              & all(coefficients_rsq_brainReg(:,idx_freq)<(mean_coef(idx_freq)+std_thresh_exclude*std_coef(idx_freq)),2);    
    Freq_data25(freq_coef).idx=idx_temp;
    Freq_data25(freq_coef).meanZS=mean(ZS_rsq(idx_temp,100:end),1);
    Freq_data25(freq_coef).Raster=ZS_rsq(idx_temp,100:end);
    Freq_data25(freq_coef).ROIs=ROI_WT_rsq(idx_temp,:);
    sum_Rois2(1) = sum_Rois2(1) + sum(idx_temp);
end

%Removing one fish
badFishNb=201808082;
for freq_coef = 1:10  
    idx_temp = Freq_data25(freq_coef).idx; 
    idx_temp = idx_temp&(idx_Fish_rsq~=badFishNb);
    Freq_data25(freq_coef).idx = idx_temp;
    Freq_data25(freq_coef).meanZS = mean(ZS_rsq(idx_temp,100:end),1);
    Freq_data25(freq_coef).Raster = ZS_rsq(idx_temp,100:end);
    Freq_data25(freq_coef).ROIs = ROI_WT_rsq(idx_temp,:);
    sum_Rois2(1) = sum_Rois2(1) + sum(idx_temp);
end
clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot

color_temp = cbrewer('div','Spectral',10);

%WB trace for all included freq
Fighandle = figure;
set(Fighandle, 'Position', [10, 10, 2400, 500]);  
for freq_coef = 1:10 
    plot(Freq_data25(freq_coef).meanZS(100:end)/max(Freq_data25(freq_coef).meanZS),'Color',color_temp(freq_coef,:));hold on;
end

Fighandle = figure;
set(Fighandle, 'Position', [10, 10, 2400, 500]);   
for freq_coef = 1:10 
    plot(Freq_data25(freq_coef).meanZS(3300:3600)/max(Freq_data25(freq_coef).meanZS),'Color',color_temp(freq_coef,:));hold on;
end

%Scatters of frequencies WB, 
Fighandle = figure;
set(Fighandle, 'Position', [10, 10, 2800, 1200]);  
for freq_coef = 1:10
    %figure;
   scatter(Freq_data25(freq_coef).ROIs(:,1),Freq_data25(freq_coef).ROIs(:,2),20,color_temp(freq_coef,:),'filled');hold on;axis([400 1400 0 600]);alpha(.5);%hold off;
end

%WB sagital all freq
Fighandle = figure;
set(Fighandle, 'Position', [10, 10, 2800, 1200]);  
for freq_coef = 1:10
   scatter(Freq_data25(freq_coef).ROIs(:,1),Freq_data25(freq_coef).ROIs(:,3),20,color_temp(freq_coef,:),'filled');hold on;axis([400 1400 0 300]);
end

ROI_all = [];
node_all = [];
for freq_coef = 1:10
   ROI_all = vertcat(ROI_all, Freq_data25(freq_coef).ROIs);
   node_all = vertcat(node_all, ones(length(Freq_data25(freq_coef).ROIs),1)*freq_coef);
end

hind_top = find(ROI_all(:,2) > 300);
[coeff,score,~,~,explained,~] = pca(ROI_all(hind_top,:));

figure;
scatter(score(:,1),score(:,2),20,node_all(hind_top),'filled'); colormap(color_temp);
%%
RegionList = {'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Stratum'};

std_thresh3=2.5;
std_thresh_exclude=2;
Freq_data_perregion2=struct();

for region_nb = 1:length(RegionList)
    region = RegionList{region_nb};
    idx_brainRegion = PerBrainRegions.(region).idx;
    idx_freq_brain = ismember(idx_brainRegion,idx_rsq);
    coefficients_rsq_brainReg = coefficients(idx_brainRegion(idx_freq_brain),:);
    
    for freq_coef = 1:10
        idx_freq = [2:1:11];
        idx_freq(freq_coef)=[];
        idx_temp = (coefficients_rsq_brainReg(:,freq_coef+1)>(mean_coef(freq_coef+1)+std_thresh3*std_coef(freq_coef+1))) ...
                   & all(coefficients_rsq_brainReg(:,idx_freq)<(mean_coef(idx_freq)+std_thresh_exclude*std_coef(idx_freq)),2);
        temp = idx_brainRegion(idx_freq_brain);
        idx_temp = temp(idx_temp);
        
        Freq_data_perregion2(freq_coef).(region).idx = idx_temp;
        Freq_data_perregion2(freq_coef).(region).meanZS = mean(ZS(idx_temp,100:3280),1);
        Freq_data_perregion2(freq_coef).(region).Raster = ZS(idx_temp,100:3280);
        Freq_data_perregion2(freq_coef).(region).ROIs = ROI_WT(idx_temp,:);
    end
    clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot
end
color_temp=cbrewer('div','Spectral',10);

%Removing one fish
for region_nb = 1:length(RegionList)    
    region = RegionList{region_nb};
    for freq_coef = 1:10       
        idx_temp = Freq_data_perregion2(freq_coef).(region).idx;
        idx_temp(idx_Fish(idx_temp)==badFishNb)=[];        
        Freq_data_perregion(freq_coef).(region).idx = idx_temp;
        Freq_data_perregion2(freq_coef).(region).meanZS = mean(ZS(idx_temp,100:3280),1);
        Freq_data_perregion2(freq_coef).(region).Raster = ZS(idx_temp,100:3280);
        Freq_data_perregion2(freq_coef).(region).ROIs = ROI_WT(idx_temp,:);
    end
    clearvars rows counter i idx_temp idx_temp2 idx_gen Fighandle yplot
end

for region_nb = 1:length(RegionList)
    region = RegionList{region_nb};
    Fighandle=figure;
    set(Fighandle, 'Position', [10, 10, 2400, 500]);
    nb_ROIs=[];
    for freq_coef=1:10
        plot(Freq_data_perregion2(freq_coef).(region).meanZS(300:end)/max(Freq_data_perregion2(freq_coef).(region).meanZS),'Color',color_temp(freq_coef,:));hold on;
        nb_ROIs=[nb_ROIs length(Freq_data_perregion(freq_coef).(region).idx)];
    end
    title(strcat(region,', nb = ',num2str(nb_ROIs)));
    print(Fighandle,strcat('..\Figures\','mean_zs_per_tone_300-3000',region,'_crit1.svg'),'-dsvg','-r0');
end

%%Region figures
for region_nb = 1:length(RegionList)
    region = RegionList{region_nb};
    Fighandle=figure;
    set(Fighandle, 'Position', [10, 10, 2400, 500]);
    for freq_coef=1:10
        scatter3(Freq_data_perregion2(freq_coef).(region).ROIs(:,1),Freq_data_perregion2(freq_coef).(region).ROIs(:,2),Freq_data_perregion2(freq_coef).(region).ROIs(:,3),40,color_temp(freq_coef,:),'filled');hold on;title(region);alpha(.8);
    end 
end
