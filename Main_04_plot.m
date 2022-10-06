
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
for k=1:19
    Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    for  i_file=1:DataInfo.condition(k,2)
        Process.loadMat=[];
        for k1=1:DataInfo.condition(k,2)
            Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
        end
        Process.DataBase=importdata(fullfile(DataInfo.location{k},char(Process.loadMat(i_file))));
        Process.DataBase=V2Pa(Process.DataBase,Process.sensorParameter.kulite_transform_ab);
        %1s数据按照整圈分成4份（截断平均）
        [Pulse,Rotor_Speed]=keyRotation_RealTime(Process.DataBase(:,end),Process.sensorParameter.fs);
        amplitude_freq_tmp=[];
        i_mean=40;
        for kk=1:i_mean
            singal=Process.DataBase(Pulse(1+floor((length(Pulse)-1)/i_mean)*(kk-1)):Pulse(floor((length(Pulse)-1)/i_mean)*kk),Process.sensorParameter.object(sensorLoc));
            sst1_tmp=resample(singal,Process.fsResample,Process.sensorParameter.fs);
            [the_freq_tmp,amplitude_freq_tmp(:,:,kk)]=fftPlot_feature(sst1_tmp,Process.fsResample,2.56,DataInfo.condition(k,1));
        end
        Process.sst1{k,i_file}=mean(sst1_tmp,3);
        the_freq{k,i_file}=the_freq_tmp;
        amplitude_freq{k,i_file}=mean(amplitude_freq_tmp,3);
        %计算频谱能量
        nlabel=find(the_freq{k,i_file}<RIband(2) & the_freq{k,i_file}>RIband(1));
        UI_FFT(k,i_file,:)=sum(amplitude_freq{k,i_file}(nlabel,:))*(the_freq{k,i_file}(2)-the_freq{k,i_file}(1))/length(nlabel);
    end
end

save NEW_RI_monitoring_DataBase_resample_100_RIband_10_22.mat


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
% 
% jet_color2=colormap(jet(15));
% h2=figure
% for k=3:15
%     plot(UI_WAVELET(k,DataInfo.condition(k,2)-80:DataInfo.condition(k,2)-80+78,2),'Color',jet_color2(k,:),'LineWidth',2)
%     hold on
% end
% scale=1;
% plot([1:67]+14,1.4*smooth(UI_WAVELET(16,1:67,2)),'Color','k','LineWidth',2)
% plot([1:76]+4, 1.4*smooth(UI_WAVELET(17,1:76,2)),'Color','b','LineWidth',2)
% plot([1:75]+5, 1.4*smooth(UI_WAVELET(18,1:75,2)),'Color','y','LineWidth',2)
% legend2 = legend('7000','7500','8000','8500','9000','9500','10000','10500','11000','11500','12000','12500','13000','14000','15000','16000')
% set(legend2,'Location','northwest');
% % close all


%导入或者存储数据
% save(['RI_monitoring_DataBase_resample_',num2str(resamplePoint),'_RIband_',num2str(RIband(1)),'_',num2str(RIband(end)),'.mat'])
% load(['RI_monitoring_DataBase_resample_',num2str(resamplePoint),'_RIband_',num2str(RIband(1)),'_',num2str(RIband(end)),'.mat'])




% close all


%% 不同转速的结果对比
h2=figure('Visible', 'on');
% set desired output size
set(h2, 'Units','centimeters')
height = 18*1.3;
width = 25*1.3;

% the last two parameters of 'Position' define the figure size
set(h2, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'InvertHardcopy', 'off',...
    'Renderer','painters'...     %recommended if there are no alphamaps
    );
set(gcf, 'color', 'white');


jet_color2=colormap(jet(15));
A=[];
for k=3:15
    A(k,:)=UI_FFT(k,DataInfo.condition(k,2)-89:DataInfo.condition(k,2)-89+78,2);
    plot(A(k,:),'Color',jet_color2(k,:),'LineWidth',2)
    hold on
end
% 
% scale=1;
% plot([1:67]+14,1.4*smooth(UI_FFT(16,1:67,2)),'Color','k','LineWidth',2)
% plot([1:76]+4, 1.4*smooth(UI_FFT(17,1:76,2)),'Color','b','LineWidth',2)
% plot([1:75]+5, 1.4*smooth(UI_FFT(18,1:75,2)),'Color','y','LineWidth',2)

legend1 = legend('7000','7500','8000','8500','9000','9500','10000','10500','11000','11500','12000','12500','13000','14000','15000','16000')
set(legend1,'Location','northwest');
famen=90:-5:29;
xuhao=round(famen*DataInfo_a(1)+DataInfo_b(1))-10;
set(gca,'FontSize',14,'XGrid','on','XTick',xuhao,...
    'XTickLabel',mat2cell(famen.',ones(length(famen),1)));
grid on
% 创建 ylabel
% ylabel({'PI（）'});
ylabel('$UI=\sum_{f1}^{f2}(W(fk)))/(f2-f1)[pa^2]$','interpreter','latex','FontSize',20)
% 创建 xlabel
xlabel('磁阀开度(%)','FontSize',30);
% 设置其余坐标区属性
title('不同转速的UI','FontSize',14)
set(gca,'FontSize',14);
saveas(h2,[save_directory,'/','FIG-overlap-50-所有转速','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')



%% 做出频带的傅里叶频谱、不同位置的UI、不同磁阀开度的UI联合图。
sensorLabel=2;
for i_rotorspeed=3:15 %固定一个转速
h1=figure('Visible', 'on');
% set desired output size
set(h1, 'Units','centimeters')
height = 12*1.3;
width = 25*1.3;

% the last two parameters of 'Position' define the figure size
set(h1, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'InvertHardcopy', 'off',...
    'Renderer','painters'...     %recommended if there are no alphamaps
    );
set(gcf, 'color', 'white');



fig=subplot('position',[0.1 0.5 0.3 0.4]);
jet_color=colormap(jet(100));
% for i_rotorspeed=3 %:size(DataInfo.condition,1) 暂时固定一个转速
    for  i_file=1:DataInfo.condition(i_rotorspeed,2)
         amplitude_freq_dB=20*log10(sqrt(amplitude_freq{i_rotorspeed,i_file}(:,sensorLabel))./2e-5);
         the_freq_temp=the_freq{i_rotorspeed,i_file};
         plot(the_freq_temp,amplitude_freq_dB,'.-','LineWidth',1,'MarkerSize',1,'Color',jet_color(i_file,:))
         hold on
    end
% end

nlabel_1BPFplot=find(the_freq_temp<30 & the_freq_temp>28);
nlabel_Transitionbandplot=find(the_freq_temp<3 & the_freq_temp>2);


colormap(gca,'jet')
ylabel('幅值/dB','FontSize',20);
% xlabel({'阶次'},'FontSize',20);
xlim([0 40])
ylim([45,120]);
txt = {['转速',num2str(DataInfo.condition(i_rotorspeed,1))],...
    ['采样',num2str(resamplePoint)],['位置',DataInfo.sensorArray{sensorLabel}]};
annotation('textbox',...
    [0.11 0.77 0.1 0.1],...
    'String',txt,...
    'FontSize',13,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[1 1 0],...
    'LineWidth',2,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
set(gca,'FontSize',15);
hold off
colorbar(gca,'TickLabels',{'95%','90%','85%','80%','75%','70%','65%','60%','55%','50%','45%','35%','30%'});
for k=1:length(RIband)
    line([RIband(k),RIband(k)],[40,120],'linestyle','--','Color','k');
end
title('不同频带的FFT幅值谱','FontSize',16)


%局部放大-且平滑处理
subplot('position',[0.1 0.11 0.3 0.3])
for  i_file=1:DataInfo.condition(i_rotorspeed,2)
    amplitude_freq_dB=20*log10(sqrt(amplitude_freq{i_rotorspeed,i_file}(:,sensorLabel))./2e-5);
    the_freq_temp=the_freq{i_rotorspeed,i_file};
    plot(the_freq_temp,smooth(amplitude_freq_dB,60),'.-','LineWidth',1,'MarkerSize',1,'Color',jet_color(i_file,:))
    hold on
end


maxRIBand1=find(the_freq_temp<RIband(2) & the_freq_temp>RIband(1));
maxRIAmplitude=max(amplitude_freq_dB(maxRIBand1));

for k=1:length(RIband)
    line([RIband(k),RIband(k)],[10,maxRIAmplitude*1.2],'linestyle','--','Color','k');
end
set(gca,'FontSize',14);
grid on
xlim([6 25])
ylim([60,maxRIAmplitude*1.1])
ylabel('局部放大(平滑处理)','FontSize',14)
xlabel('阶次','FontSize',14)
hold off

%不同位置的UI(频带内的能量)
subplot('position',[0.46 0.11 0.15 0.8])
jet_color=colormap(jet(100));
for  i_file=1:DataInfo.condition(i_rotorspeed,2)
    Level=1000;
    plot(Level*permute(UI_FFT(i_rotorspeed,i_file,:),[3,1,2]),DataInfo.sensorLoc,'+-','LineWidth',2,'Color',jet_color(i_file,:))
    hold on
end
% famen=[30:1:75];
% xuhao=round((famen-DataInfo.b(label))./DataInfo.a(label));
% maximum1=max(max(PI1));
% jet_color1=colormap(jet(length(famen)));
% for k=1:length(famen)
%     plot(PI1(xuhao(k),:)/maximum1,DataInfo.sensorLoc,'+-','LineWidth',2,'Color',jet_color1(length(famen)-k+1,:));
%     hold on
% end
grid on
set(gca,'FontSize',14,'YGrid','on','YTick',DataInfo.sensorLoc,...
    'YTickLabel',DataInfo.sensorArray);
set(gca,'FontSize',14);
xlabel('单位频带能量/pa^2x10^-^3','FontSize',14);
title('不同位置的UI','FontSize',16)
% xlim([0 1])
ylim([-1 101]);

%
subplot('position',[0.67 0.11 0.31 0.8])

plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),10),'LineWidth',2,'Color','k')
hold on
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),9),'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),8),'LineWidth',2,'Color',[0.4940 0.1840 0.5560])
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),7),'LineWidth',2,'Color',[0.3010 0.7450 0.9330])
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),6),'LineWidth',2,'Color',[96 96 96]/255)
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),5),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),4),'LineWidth',2,'Color','y')
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),3),'LineWidth',2,'Color','b')
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),2),'LineWidth',2,'Color','g')
plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),1),'LineWidth',2,'Color','r')
hold off

maxUI_FFT=max(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),2));
%添加SI虚线
line([DataInfo.condition(i_rotorspeed,2),DataInfo.condition(i_rotorspeed,2)],[0,maxUI_FFT*1.2],'linestyle','--','Color','k');
ylim([[0,maxUI_FFT*1.05]]);
legend1 = legend('C1','R8','R7','R6','R5','R4','R3','R2','R1','B1','SI')
set(legend1,'NumColumns',1,'Location','northwest');

%利用序号和阀门开度的关系换算(xuhao=(famen-b)/a;)
famen=90:-5:29;
xuhao=round(famen*DataInfo_a(i_rotorspeed)+DataInfo_b(i_rotorspeed));
set(gca,'FontSize',14,'XGrid','on','XTick',xuhao,...
    'XTickLabel',mat2cell(famen.',ones(length(famen),1)));
grid on
% 创建 ylabel
% ylabel({'PI（）'});
ylabel('$UI=\sum_{f1}^{f2}(W(fk)))/(f2-f1)[pa^2]$','interpreter','latex','FontSize',20)
% 创建 xlabel
xlabel('磁阀开度(%)','FontSize',30);
% 设置其余坐标区属性
title('不同磁阀开度的UI','FontSize',14)
set(gca,'FontSize',14);
saveas(h1,[save_directory,'/','FIG-overlap-50-转速',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')
end









%%



% saveas(h1,[save_directory,'/','FIG-转速',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'.png'])
% saveas(h1,[save_directory,'/','FIG-转速',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'.fig'])
% saveas(h1,[save_directory,'/','FIG-转速',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'.eps'])

% cleanfigure
% matlab2tikz([save_directory,'/','FIG-转速',num2str(RotorSpeed),'rpm','-采样率',num2str(resamplePoint),'.tex']);









