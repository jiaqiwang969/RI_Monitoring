%出图代码，Figure：傅里叶和小波变换结果
%12000rpm
%三种不同的方法对比
%采样率116点/圈
%RIband=[8:24]
%仅对比B1和R1

clc
clear
close all


%% 参数
for N_caiyang=3:6;%每个叶道的点数
    resamplePoint=N_caiyang*29; %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
    RIband=[8:24];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
    sensorLoc=[1:2];




    %%
    addpath(genpath('subfunction'));
    addpath(genpath('Test'));
    addpath(genpath('synchrosqueezing'));



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
DataInfo.condition(3,:) =[7000;95-1;25;95];
DataInfo.condition(4,:) =[7500;93-1;25+1/3;95];
DataInfo.condition(5,:) =[8000;93-1;25+1/3*2;95];
DataInfo.condition(6,:) =[8500;95-1;25+1/3*3;95];
DataInfo.condition(7,:) =[9000;95-1;25+1/3*4;95];
DataInfo.condition(8,:) =[9500;95-1;25+1/3*5;95];
DataInfo.condition(9,:) =[10000;95-1;25+1/3*6;95];
DataInfo.condition(10,:)=[10500;93-1;25+1/3*7;95];
DataInfo.condition(11,:)=[11000;94-1;25+1/3*8;95];
DataInfo.condition(12,:)=[11500;97-1;25+1/3*9;95];
DataInfo.condition(13,:)=[12000;93-1;25+1/3*10;95];
DataInfo.condition(14,:)=[12500;91-1;25+1/3*11;95];
DataInfo.condition(15,:)=[13000;92-1;29;95];
DataInfo.condition(16,:)=[14000;67-1;29;70];
DataInfo.condition(17,:)=[15000;76-1;29;70];
DataInfo.condition(18,:)=[16000;75-1;29;70];




    %目的：统一选择磁阀开度范围（横坐标固定）
    %计算磁阀开度和序号的线性对应关系(famen=b-xuhao*a;)
    x1=DataInfo.condition(:,4);x2=DataInfo.condition(:,3);
    y1=1;y2=DataInfo.condition(:,2);
    DataInfo_a=(y2-y1)./(x2-x1);
    DataInfo_b=y1-DataInfo_a.*x1;

    % 完善工作 利用上述数据信息，构造一个与UI_FFT同等长度的矩阵18*96
    XUHAO=repmat([1:max(DataInfo.condition(:,2))],size(DataInfo.condition,1),1);
    FAMEN=(XUHAO-DataInfo_b)./DataInfo_a;

    %%

    for k=1:length(DataInfo.condition)
        DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(DataInfo.condition(k,1))];
        DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
    end

    DataInfo.sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
    %10个传感器的相对弦长位置
    DataInfo.sensorLoc=[-1;5;13;23;35;47;58;71;82;101];


    %转速的影响-FFT
    the_freq=[];
    amplitude_freq=[];
    UI_FFT=[];
    Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));


    for k=13
        Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
        Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
        parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
            % %         Process.loadMat=[];
            %         for k1=1:DataInfo.condition(k,2)
            %             Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
            %         end
            DataBase=importdata(fullfile(DataInfo.location{k},['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(i_file),'.mat']));
            DataBase=V2Pa(DataBase,Process.sensorParameter.kulite_transform_ab);

            % 执行等角度采样操作:完全按照阶次谱来执行
            % n_round_Point=(Pulse(n_round+1)-Pulse(n_round))/29;
            % xuhao=Pulse(n_round)+(n_blade-1)*n_round_Point+(1/(2*N_caiyang)+n_caiyang/N_caiyang)*n_round_Point
            QUAN=112; %圈数按照最小转速(7000rpm)统一固定为112，目的是为了后续统一频率分布
            [Pulse,Rotor_Speed]=keyRotationRealTime(DataBase(:,end),Process.sensorParameter.fs);
            xuhao=[];
            for n_round=1:QUAN
                n_round_Point=(Pulse(n_round+1)-Pulse(n_round))/29;
                for n_blade=1:29
                    for n_caiyang=1:N_caiyang
                        xuhao=[xuhao round(Pulse(n_round)+(n_blade-1)*n_round_Point+(1/(2*N_caiyang)+n_caiyang/N_caiyang)*n_round_Point)];
                    end
                end
            end
            sst1All{N_caiyang,i_file}=DataBase(xuhao,Process.sensorParameter.object(sensorLoc));
            %%利用流动失稳的循环平稳特性分离出tonal noise和broadband noise
            %背景流
            data_tonal_rms=permute(mean(reshape(sst1All{N_caiyang,i_file},QUAN,29*N_caiyang,length(Process.sensorParameter.object(sensorLoc))),1),[2,3,1]);
            data_tonal{N_caiyang,i_file}=kron(ones(QUAN,1),data_tonal_rms);
            data_broadband{N_caiyang,i_file}=sst1All{N_caiyang,i_file}-data_tonal{N_caiyang,i_file};
            %尝试减去方差！！看UI曲线的变化
            %         data_broadband{N_caiyang,i_file}= data_broadband{N_caiyang,i_file}./std(data_broadband{N_caiyang,i_file});

            %       data_diff{N_caiyang,i_file}=reshape(diff(reshape(sst1All{N_caiyang,i_file},QUAN,N_caiyang*29,length(Process.sensorParameter.object(sensorLoc))),1),(QUAN-1)*N_caiyang*29,length(Process.sensorParameter.object(sensorLoc)));

            %% 方法1:不平均
            P2 = 2*abs(fft(data_broadband{N_caiyang,i_file})/length(data_broadband{N_caiyang,i_file}));
            the_freq1{N_caiyang,i_file}= N_caiyang*29*(0:(length(data_broadband{N_caiyang,i_file})/2.56))/length(data_broadband{N_caiyang,i_file});
            amplitude_freq1{N_caiyang,i_file}= P2(1:length(data_broadband{N_caiyang,i_file})/2.56+1,:);
            nlabel1=find(the_freq1{N_caiyang,i_file}<RIband(end) & the_freq1{N_caiyang,i_file}>RIband(1));
            UI_FFT_noAverage(N_caiyang,i_file,:)=sum(amplitude_freq1{N_caiyang,i_file}(nlabel1,:).^2)*(the_freq1{N_caiyang,i_file}(2)-the_freq1{N_caiyang,i_file}(1));

        end
        %
        %     %% 方法2:stft：平均+overlap
        %     %
        %         for iSensor=1:length(Process.sensorParameter.object(sensorLoc))
        %             parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
        %                 %various options and parameters
        %                 CWTopt=struct('gamma',eps,'type','morlet','mu',6,'s',2,'om',0,'nv',64,'freqscale','linear');
        %                 STFTopt=struct('gamma',eps,'type','gauss','mu',1,'s',2.5,'om',0.5,'winlen',256*4,'squeezing','full');
        %                 %Short-time Fourier transform (STFT)
        %                 dt=1/(N_caiyang*29);
        %                 [Sx,fs,dSx] = stft_fw(data_broadband{N_caiyang,i_file}(:,iSensor), dt, STFTopt);
        %                 the_freq2{N_caiyang,i_file}= 2*[0:size(Sx,1)-1]/size(Sx,1)*29;
        %                 amplitude_freq2x{N_caiyang,i_file}=mean(abs(Sx)/2,2);  %这里除以2，才使得PI结果统一？
        %
        %                 %         figure
        %                 %         plot(the_freq1{N_caiyang,i_file},amplitude_freq1{N_caiyang,i_file}(:,2));
        %                 %         hold on
        %                 %         plot(the_freq2{N_caiyang,i_file},amplitude_freq2{N_caiyang,i_file});
        %                 %
        %                 %         % 用于验证信号
        %                 %        xNew = stft_iw(Sx, fs, STFTopt).';
        %                 %        figure;tplot(Sx, 1:length(data_broadband{N_caiyang,i_file}(:,2)), fs); colorbar; title('STFT','FontSize',20); xlabel('Time (seconds)','FontSize',20); ylabel('Frequency (hz)', 'FontSize',20);
        %                 %        figure(); plot(1:length(data_broadband{N_caiyang,i_file}(:,2)),[data_broadband{N_caiyang,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
        %
        %                 % % STFT Synchrosqueezing transform
        %                 %
        %                 % [Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(1:length(data_broadband{N_caiyang,i_file}(:,2)), data_broadband{N_caiyang,i_file}(:,2), STFTopt);
        %                 % xNew = synsq_stft_iw(Tx, fs, STFTopt).';
        %                 % figure(); tplot(Tx, 1:length(data_broadband{N_caiyang,i_file}(:,2)), fs); colorbar; title('STFT Synchrosqueezing','FontSize',20); xlabel('Time (seconds)','FontSize',20); ylabel('Frequency (hz)', 'FontSize',20);
        %                 % figure(); plot(1:length(data_broadband{N_caiyang,i_file}(:,2)),[data_broadband{N_caiyang,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
        %
        %
        %                 %%
        %                 nlabel2=find(the_freq2{N_caiyang,i_file}<RIband(end) & the_freq2{N_caiyang,i_file}>RIband(1));
        %                 UI_FFT_stftAverage(k,i_file,iSensor)=sum(amplitude_freq2x{N_caiyang,i_file}(nlabel2).^2)*(the_freq2{N_caiyang,i_file}(2)-the_freq2{N_caiyang,i_file}(1));
        %             end
        %             for i_file=1:DataInfo.condition(k,2)
        %                 amplitude_freq2{N_caiyang,i_file}(:,iSensor)=amplitude_freq2x{N_caiyang,i_file};
        %             end
        %         end
        %
        % RIband=[8:22];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
        %
        %
        %     %% 方法3:wavelet：平均+overlap
        %
        %     for iSensor=1:length(Process.sensorParameter.object(sensorLoc))
        %         parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
        %             fk=[1:0.5:58];
        %             pad = 0;                                   % 2次幂0填充
        %             dj = 0.125;                                % this will do 4 sub-octaves per octave
        %             j1 = 200;                                  % this says do 7 powers-of-two with dj sub-octaves each
        %             lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
        %             mother = 'morlet';
        %             fsResample=DataInfo.condition(k,1)/60*N_caiyang*29;
        %             dt = 1/fsResample ;
        %             s0 = 15*dt;
        %
        %             period=1./fk/(DataInfo.condition(k,1)/60);
        %             % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
        %             fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
        %             scale=period/fourier_factor;
        %             fk=1./period/(DataInfo.condition(k,1)/60);
        %
        %             [wave,period] = wavelet_feature(data_broadband{N_caiyang,i_file},dt,pad,dj,mother,scale);
        %             the_freq3{N_caiyang,i_file}= fk;
        %             amplitude_freq3{N_caiyang,i_file}= permute(mean(abs(wave),2),[1,3,2]);
        %             nlabel3=find(the_freq3{N_caiyang,i_file}<RIband(end) & the_freq3{N_caiyang,i_file}>RIband(1));
        %             UI_FFT_waveletAverage(k,i_file,:)=sum(amplitude_freq3{N_caiyang,i_file}(nlabel3,:).^2)*(the_freq3{N_caiyang,i_file}(2)-the_freq3{N_caiyang,i_file}(1));
        %         end
        %     end
        % %
        % end

    end
end



%% 做出频带的傅里叶频谱、不同位置的UI、不同磁阀开度的UI联合图。
    RIband1=[1:80]
    RIbandName=[];
    RIbandName{10}=10;
    RIbandName{20}=20;
    RIbandName{30}=30;
    RIbandName{40}=40;


amplitude_freq=amplitude_freq1;
UI_FFT=UI_FFT_noAverage;
the_freq=the_freq1;

i_rotorspeed=13;
h1=figure('Visible', 'on');
% set desired output size
set(h1, 'Units','centimeters')
height = 16*1.3;
width = 25*1.3;

% the last two parameters of 'Position' define the figure size
set(h1, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'InvertHardcopy', 'off',...
    'Renderer','painters'...     %recommended if there are no alphamaps
    );
set(gcf, 'color', 'white');



subplot('position',[0.1 0.1 0.8 0.8])

sensorLabel=2;
% for N_caiyang=1:4 %固定一个转速
    RIband=[8:24]


    %     %不同位置的UI(频带内的能量)
    %     subplot('position',[0.1 0.11 0.15 0.8])
    %     jet_color=colormap(jet(100));
    %     for  i_file=1:DataInfo.condition(i_rotorspeed,2)
    %         Level=1/1000;
    %         plot(Level*permute(UI_FFT(i_rotorspeed,i_file,:),[3,1,2]),DataInfo.sensorLoc,'+-','LineWidth',2,'Color',jet_color(i_file,:))
    %         hold on
    %     end
    %     % famen=[30:1:75];
    %     % xuhao=round((famen-DataInfo.b(label))./DataInfo.a(label));
    %     % maximum1=max(max(PI1));
    %     % jet_color1=colormap(jet(length(famen)));
    %     % for k=1:length(famen)
    %     %     plot(PI1(xuhao(k),:)/maximum1,DataInfo.sensorLoc,'+-','LineWidth',2,'Color',jet_color1(length(famen)-k+1,:));
    %     %     hold on
    %     % end
    %     grid on
    %     set(gca,'FontSize',20,'YGrid','on','YTick',DataInfo.sensorLoc,...
    %         'YTickLabel',DataInfo.sensorArray);
    %     set(gca,'FontSize',20);
    %     xlabel('频带能量/pa^2x10^3','FontSize',20);
    %     title('不同位置的I_R_I','FontSize',16)
    %     % xlim([0 1])
    %     ylim([-1 101]);

    %
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),10),'LineWidth',2,'Color','k')
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),9),'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),8),'LineWidth',2,'Color',[0.4940 0.1840 0.5560])
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),7),'LineWidth',2,'Color',[0.3010 0.7450 0.9330])
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),6),'LineWidth',2,'Color',[96 96 96]/255)
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),5),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),4),'LineWidth',2,'Color','y')
    %     plot(UI_FFT(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),3),'LineWidth',2,'Color','b')
%     plot(UI_FFT_noAverage(1,1:DataInfo.condition(i_rotorspeed,2),2),'*-','LineWidth',2)

%     plot(UI_FFT_noAverage(2,1:DataInfo.condition(i_rotorspeed,2),2),'+-','LineWidth',2)
    plot(UI_FFT_noAverage(3,1:DataInfo.condition(i_rotorspeed,2),2),'^-','LineWidth',2)
        hold on

    plot(UI_FFT_noAverage(4,1:DataInfo.condition(i_rotorspeed,2),2),'o-','LineWidth',2)
    plot(UI_FFT_noAverage(5,1:DataInfo.condition(i_rotorspeed,2),2),'+-','LineWidth',2)
    plot(UI_FFT_noAverage(6,1:DataInfo.condition(i_rotorspeed,2),2),'*-','LineWidth',2)
%     plot(UI_FFT_noAverage(7,1:DataInfo.condition(i_rotorspeed,2),2),'+-','LineWidth',2)
    txt = {...
        ['转速',num2str(DataInfo.condition(i_rotorspeed,1))], ['位置',DataInfo.sensorArray{2}],['RI带',num2str(RIband(1)),'-',num2str(RIband(end))]};
    annotation('textbox',...
        [0.2 0.77 0.17 0.13],...
        'String',txt,...
        'FontSize',20,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);

%     plot(UI_FFT_noAverage(N_caiyang,1:DataInfo.condition(i_rotorspeed,2),1),'-','LineWidth',2,'Color','r')
%     txt = {['转速',num2str(DataInfo.condition(i_rotorspeed,1))],...
%         ['采样',num2str(N_caiyang*29)],...
%         ['RI带',num2str(RIband(1)),'-',num2str(RIband(end))]};
    %     annotation('textbox',...
    %         [0.11+0.1 0.76 0.1 0.13],...
    %         'String',txt,...
    %         'FontSize',13,...
    %         'FontName','Arial',...
    %         'LineStyle','--',...
    %         'EdgeColor',[1 1 0],...
    %         'LineWidth',2,...
    %         'BackgroundColor',[0.9  0.9 0.9],...
    %         'Color',[0.84 0.16 0]);
    %     set(gca,'FontSize',15);
    %     plot(1/2400*UI_FFT_waveletAverage(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),2),'+-','LineWidth',2,'Color','g')
    %     plot(1/2400*UI_FFT_waveletAverage(i_rotorspeed,1:DataInfo.condition(i_rotorspeed,2),1),'+-','LineWidth',2,'Color','r')




    maxUI_FFT=max(UI_FFT_noAverage(N_caiyang,1:DataInfo.condition(i_rotorspeed,2),2));
    %添加SI虚线
    %     line([DataInfo.condition(i_rotorspeed,2),DataInfo.condition(i_rotorspeed,2)],[0,maxUI_FFT*1.2],'linestyle','--','Color','k');
    %     ylim([[0,maxUI_FFT*1.05]]);



    %     saveas(h1,[save_directory,'/','FIG-无平均-方差',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')

% end


%利用序号和阀门开度的关系换算(xuhao=(famen-b)/a;)
famen=90:-5:29;
xuhao=round(famen*DataInfo_a(i_rotorspeed)+DataInfo_b(i_rotorspeed));
set(gca,'FontSize',20,'XGrid','on','XTick',xuhao,...
    'XTickLabel',mat2cell(famen.',ones(length(famen),1)));
grid on
% 创建 ylabel
% ylabel({'PI（）'});
ylabel('FFT频带能量E_R_I','FontSize',20)
% 创建 xlabel
xlabel('磁阀开度(%)','FontSize',20);
% 设置其余坐标区属性
% title('','FontSize',20)
set(gca,'FontSize',20);
    legend1 = legend('87','116','145','174')
    set(legend1,'NumColumns',1,'Location','northwest','FontSize',20);
         box off
    param = struct('XLabel', '磁阀开度(%)', 'YLabel', '$pa^2$');
    DrawAxisWithArrow(gca, param);

% subplot('position',[0.55 0.11 0.39 0.8])
% 
% for N_caiyang=1:7;%每个叶道的点数
%     plot(the_freq1{N_caiyang,73},200*N_caiyang+100*log10(amplitude_freq1{N_caiyang,73}(:,2)/2e-5)-130);
%     hold on
% end
% ylim([500 2200])
% grid on
% xlabel('阶次','FontSize',30);
% 
% 
% % for k=1:length(RIband)
% %         line([RIband(k),RIband(k)],[0,2200],'linestyle','--','Color','k');
% % end
% 
%      line([8,8],[0,21050],'linestyle','--','Color','k');
%      line([25,25],[0,2150],'linestyle','--','Color','k');
%       text(5,2160,'RI带下限')
%       text(20,2160,'RI带上限')
%     set(gca,'FontSize',15);
%         set(gca,'XGrid','on','XTick',RIband1,...
%         'XTickLabel',RIbandName);
% 
% title('不同采样率的频谱图','FontSize',20)
% 
% 
% set(gca,'FontSize',20,'YGrid','on','YTick',[700 900 1100 1300 1500 1700 1900],...
%     'YTickLabel',{'29','58','87','116','145','174','203'});

saveas(h1,[save_directory,'/','FIG04-caiyanglv-',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'eps')
saveas(h1,[save_directory,'/','FIG04-caiyanglv-',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'fig')
saveas(h1,[save_directory,'/','FIG04-caiyanglv-',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')

















