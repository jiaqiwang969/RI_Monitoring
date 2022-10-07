
clc
clear
close
addpath(genpath('subfunction'));
%% 保存图像至指定文件夹

save_directory = ['Result-',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end

%DataInfo:包含数据的位置、数据转速、数据序列、对应阀门开度、传感器位置等信息
%RotorSpeed、xuhao_end;famen_init、famen_end
%RotorSpeed、xuhao_end;famen_init、famen_end
DataInfo.condition(1,:) =[6000;92-1;29;95];
DataInfo.condition(2,:) =[6500;94-1;29;95];
DataInfo.condition(3,:) =[7000;95-1;29;95];
DataInfo.condition(4,:) =[7500;93-1;29;95];
DataInfo.condition(5,:) =[8000;93-1;29;95];
DataInfo.condition(6,:) =[8500;95-1;29;95];
DataInfo.condition(7,:) =[9000;95-1;29;95];
DataInfo.condition(8,:) =[9500;95-1;29;95];
DataInfo.condition(9,:) =[10000;95-1;29;95];
DataInfo.condition(10,:)=[10500;93-1;29;95];
DataInfo.condition(11,:)=[11000;94-1;29;95];
DataInfo.condition(12,:)=[11500;97-1;29;95];
DataInfo.condition(13,:)=[12000;93-1;29;95];
DataInfo.condition(14,:)=[12500;91-1;29;95];
DataInfo.condition(15,:)=[13000;92-1;29;95];
DataInfo.condition(16,:)=[14000;67-1;29;67];
DataInfo.condition(17,:)=[15000;76-1;29;67];
DataInfo.condition(18,:)=[16000;75-1;29;67];

%目的：统一选择阀门开度范围（横坐标固定）
%计算阀门开度和序号的线性对应关系(famen=xuhao*a+b;)
DataInfo.a=(DataInfo.condition(:,4)-DataInfo.condition(:,3))./(DataInfo.condition(:,2)-1);
DataInfo.b=DataInfo.condition(:,3)-1.*DataInfo.a;

for k=1:length(DataInfo.condition)
    DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(DataInfo.condition(k,1))];
    DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
end

DataInfo.sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
%10个传感器的相对弦长位置
DataInfo.sensorLoc=[-1;5;13;23;35;47;58;71;82;101];
DataInfo.label=75:-5:29;
DataInfo.xuhao=round((DataInfo.label-DataInfo.b)./DataInfo.a);
for k=1:length(DataInfo.label)
    DataInfo.labelName{k}=[num2str(DataInfo.label(k))];
end

%转速的影响-FFT
resamplePoint=100;
RIband=[9:11];
sensorLoc=[1:10];
the_freq=[];
amplitude_freq=[];
UI_FFT=[];
Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));
for k=1:18
    Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    for  i_file=1:DataInfo.condition(k,2)
        Process.loadMat=[];
        for k1=1:DataInfo.condition(k,2)
            Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
        end
        Process.DataBase=importdata(fullfile(DataInfo.location{k},char(Process.loadMat(i_file))));
        Process.DataBase=V2Pa(Process.DataBase,Process.sensorParameter.kulite_transform_ab);
        Process.sst1{k,i_file}=resample(Process.DataBase(:,Process.sensorParameter.object(sensorLoc)),Process.fsResample,Process.sensorParameter.fs);
        [Pulse,Rotor_Speed]=keyRotation_RealTime(Process.DataBase(:,end),Process.sensorParameter.fs);
        for kk=1:4 %分成4份截断平均
            [the_freq_tmp,amplitude_freq_tmp(:,:,kk)]=fftPlot_feature(Process.sst1{k,i_file}(1+floor((length(Pulse)-1)/4)*(k-1):floor((length(Pulse)-1)/4)*k,:),Process.fsResample,2.56,DataInfo.condition(k,1));
        end
        the_freq{k,i_file}=the_freq_tmp;
        amplitude_freq{k,i_file}=mean(amplitude_freq_tmp,3)
        %计算频谱能量
        nlabel=find(the_freq{k,i_file}<RIband(2) & the_freq{k,i_file}>RIband(1));
        UI_FFT(k,i_file,:)=sum(amplitude_freq{k,i_file}(nlabel,:))*(the_freq{k,i_file}(2)-the_freq{k,i_file}(1))/length(nlabel);
    end
end

jet_color2=colormap(jet(15));
h1=figure
A=[];
for k=3:15
    A(k,:)=UI_FFT(k,DataInfo.condition(k,2)-90:DataInfo.condition(k,2)-90+78,2);
    plot(smooth(A(k,:)),'Color',jet_color2(k,:),'LineWidth',2)
    hold on
end
% 
scale=1;
plot([1:67]+14,1.4*smooth(UI_FFT(16,1:67,2)),'Color','k','LineWidth',2)
plot([1:76]+4, 1.4*smooth(UI_FFT(17,1:76,2)),'Color','b','LineWidth',2)
plot([1:75]+5, 1.4*smooth(UI_FFT(18,1:75,2)),'Color','y','LineWidth',2)

legend1 = legend('7000','7500','8000','8500','9000','9500','10000','10500','11000','11500','12000','12500','13000','14000','15000','16000')
set(legend1,'Location','northwest');

%转速的影响-Wavelet
resamplePoint=100;
RIband=[13:15];
sensorLoc=[1:10];

Wavlet.fk=RIband;

Wavlet.pad = 0;                                   % 2次幂0填充
Wavlet.dj = 0.125;                                % this will do 4 sub-octaves per octave
Wavlet.j1 = 200;                                  % this says do 7 powers-of-two with dj sub-octaves each
Wavlet.lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
Wavlet.mother = 'morlet';


UI_WAVELET=[];
Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));
for k=1:18
    Process.fsResample=round(resamplePoint*DataInfo.condition(k,1)/60/10)*10;
    Process.dt = 1/Process.fsResample ;
    Wavlet.s0 = 15*Process.dt;

    Wavlet.period=1./Wavlet.fk/(DataInfo.condition(k,1)/60);
    % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
    Wavlet.fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
    Wavlet.scale=Wavlet.period/Wavlet.fourier_factor;
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    
    for  i_file=1:DataInfo.condition(k,2)
        Process.loadMat=[];
        for k1=1:DataInfo.condition(k,2)
            Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
        end
        Process.DataBase=importdata(fullfile(DataInfo.location{k},char(Process.loadMat(i_file))));
        Process.DataBase=V2Pa(Process.DataBase,Process.sensorParameter.kulite_transform_ab);
        Process.sst1{k,i_file}=resample(Process.DataBase(:,Process.sensorParameter.object(sensorLoc)),Process.fsResample,Process.sensorParameter.fs);
        %[the_freq{k,i_file},amplitude_freq{k,i_file}]=fftPlot_dB_universal_feature(Process.sst1{k,i_file},Process.fsResample,2.56,DataInfo.condition(k));

        %小波计算
        [Wavlet.wave,Wavlet.period] = wavelet_feature(Process.sst1{k,i_file},Process.dt,Wavlet.pad,Wavlet.dj,Wavlet.mother,Wavlet.scale);

        Wavlet.power = (abs(Wavlet.wave)).^2/length(Wavlet.wave)^2;        % compute wavelet power spectrum
        Wavlet.fk=1./Wavlet.period/(DataInfo.condition(k,1)/60);
        Wavlet.global_ws{k,i_file}=permute(sum(Wavlet.power,2),[1 3 2]);

        
        UI_WAVELET(k,i_file,:)=sum(Wavlet.global_ws{k,i_file})/length(RIband);
    end
end

jet_color2=colormap(jet(15));
h2=figure
for k=3:15
    plot(UI_WAVELET(k,DataInfo.condition(k,2)-80:DataInfo.condition(k,2)-80+78,2),'Color',jet_color2(k,:),'LineWidth',2)
    hold on
end
scale=1;
plot([1:67]+14,1.4*smooth(UI_WAVELET(16,1:67,2)),'Color','k','LineWidth',2)
plot([1:76]+4, 1.4*smooth(UI_WAVELET(17,1:76,2)),'Color','b','LineWidth',2)
plot([1:75]+5, 1.4*smooth(UI_WAVELET(18,1:75,2)),'Color','y','LineWidth',2)
legend2 = legend('7000','7500','8000','8500','9000','9500','10000','10500','11000','11500','12000','12500','13000','14000','15000','16000')
set(legend2,'Location','northwest');
% close all

%存储结果；存储mat文件

save RI_monitoring_DataBase_resample_100_RIband_9_11.mat
