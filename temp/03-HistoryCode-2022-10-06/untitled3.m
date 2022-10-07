
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



%DataInfo:包含数据的位置、数据转速、数据序列、对应磁阀开度、传感器位置等信息

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


%目的：统一选择磁阀开度范围（横坐标固定）
%计算磁阀开度和序号的线性对应关系(famen=b-xuhao*a;)
x1=DataInfo.condition(:,4);x2=DataInfo.condition(:,3);
y1=1;y2=DataInfo.condition(:,2);
DataInfo_a=(y2-y1)./(x2-x1);
DataInfo_b=y1-DataInfo_a.*x1;

for k=1:length(DataInfo.condition)
    DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(DataInfo.condition(k,1))];
    DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
end

DataInfo.sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
%10个传感器的相对弦长位置
DataInfo.sensorLoc=[-1;5;13;23;35;47;58;71;82;101];



%转速的影响-FFT
resamplePoint=100; %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
RIband=[11:22];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
sensorLoc=[1:10];  
the_freq=[];
amplitude_freq=[];
UI_FFT=[];
Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));
for k=13
    Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    for  i_file=67
        Process.loadMat=[];
        for k1=1:DataInfo.condition(k,2)
            Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
        end
        Process.DataBase=importdata(fullfile(DataInfo.location{k},char(Process.loadMat(i_file))));
        Process.DataBase=V2Pa(Process.DataBase,Process.sensorParameter.kulite_transform_ab);
        Process.sst1{k,i_file}=resample(Process.DataBase(:,Process.sensorParameter.object(sensorLoc)),Process.fsResample,Process.sensorParameter.fs);
        [the_freq{k,i_file},amplitude_freq{k,i_file}]=fftPlot_feature(Process.sst1{k,i_file},Process.fsResample,2.56,DataInfo.condition(k,1));
        
        %计算频谱能量
        nlabel=find(the_freq{k,i_file}<RIband(2) & the_freq{k,i_file}>RIband(1));
        UI_FFT(k,i_file,:)=sum(amplitude_freq{k,i_file}(nlabel,:))*(the_freq{k,i_file}(2)-the_freq{k,i_file}(1))/length(nlabel);
    end
end



figure
plot(  Process.DataBase(:,2))

% [the_freq,amplitude_freq]=fftPlot_feature(Process.DataBase(:,2),Process.fsResample,2.56,12000)

figure

% plot(the_freq,amplitude_freq)

figure

plot(abs(fft(Process.DataBase(:,2)))*2)








