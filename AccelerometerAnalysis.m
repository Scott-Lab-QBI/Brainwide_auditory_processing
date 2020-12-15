figure;plot(mean(AlignedSignals_x,2));
figure;plot(mean(AlignedSignals_x(25000:35000,,),2));

timing=([25000,35000;35000,45000;48000,58000;58000,70000;70000,82000;83000,93000;94000,104000;104000,114000;114000,124000;125000,135000;32000,38000]);

figure;plot(temp(32000:38000));


temp=mean(AlignedSignals_x,2);
figure;
for i=1:size(timing,1)
    plot(temp(timing(i,1):timing(i,2)));
    pause;
end

fs=3938.4615;

test_filt=temp(timing(1,1):timing(1,2));
figure;bandpass(test_filt,[80 120],fs);


figure;periodogram(test_filt,[],length(test_filt),fs)


[pxx,f] = pspectrum(test_filt,fs,'FrequencyResolution',10);

figure;plot(f,pow2db(pxx));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')


figure;bandpass(test_filt,[90 4000],fs);

[pxx,f] = pspectrum(temp(timing(2,1):timing(2,2)),fs,'FrequencyResolution',10);

figure;plot(f,pow2db(pxx));
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')


temp_ts=bandpass(test_filt,[80 4000],fs);
g=9.80665;calibration=98.7;


AlignedSignals_3D_mean=horzcat(mean(AlignedSignals_x,2),mean(AlignedSignals_y,2),mean(AlignedSignals_z,2));
AlignedSignals_3D=cat(3,AlignedSignals_x,AlignedSignals_y(:,1:5),AlignedSignals_z(:,1:5));




acceleration=zeros(10,3,5);
for ik=1:size(acceleration,2)
    for kj=1:size(acceleration,3)
        test_filt=AlignedSignals_3D(:,kj,ik);
        for ij=1:size(acceleration,1)            
            temp_ts=highpass(test_filt(timing(ij,1):timing(ij,2)),80,fs);
            acceleration(ij,ik,kj)=max(temp_ts);
        end
    end
end
acceleration=acceleration*g*calibration/100;

figure;plot(mean(acceleration,3));

PowerDB=zeros(length(pxx),size(timing,1),3,5);
for ik=1:size(PowerDB,3)
    for kj=1:size(PowerDB,4)
        test_filt=AlignedSignals_3D(:,kj,ik);
        for ij=1:size(PowerDB,2)            
            [pxx,f] = pspectrum(test_filt(timing(ij,1):timing(ij,2)),fs,'FrequencyResolution',10);
            PowerDB(:,ij,ik,kj)=pow2db(pxx);
        end
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=3;
yplot=1;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
for i=1:3
    axes(ha(i));
    imagesc(squeeze(mean(PowerDB(:,:,i,:),4)));
end
colormap hot;


Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=3;
yplot=1;
ha = tight_subplot(yplot,xplot,[.01 .01],[.01 .01],[.01 .01]);
for i=1:3
    axes(ha(i));
    plot(squeeze(mean(PowerDB(:,:,i,:),4)));
end
colormap hot;


Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1200, 1200]);
xplot=4;
yplot=3;
ha = tight_subplot(yplot,xplot,[.05 .05],[.05 .05],[.05 .05]);
for i=1:size(timing,1)
    axes(ha(i));
    plot(f,squeeze(mean(PowerDB(:,i,:,:),4)),'LineWidth',3);ylim([-100 -20]);
end
print(Fighandle,'FrequencyPowerDistribution','-dsvg');

PrismTemp=nan(size(acceleration,1),size(acceleration,2)*size(acceleration,3));
for axis_nb=1:3    
    PrismTemp(:,1+5*(axis_nb-1):5+5*(axis_nb-1))=acceleration(:,axis_nb,:);
end

