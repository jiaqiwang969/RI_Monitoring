%% 重新按照新的思路整理代码：
%暂时不做了！还是按照原来的思路来
%% 4页纸

%% 1. 初始数据
%7000:500:13000rpm
SI_position(1,:) =[7000;95;0.4925];
SI_position(2,:) =[7500;93;0.633931];
SI_position(3,:) =[8000;94;0.999];
SI_position(4,:) =[8500;95;0.438];
SI_position(5,:) =[9000;95;0.11];
SI_position(6,:) =[9500;95;0.45];
SI_position(7,:) =[10000;95;0.827764];
SI_position(8,:) =[10500;93;0.8962];
SI_position(9,:) =[11000;94;0.417];
SI_position(10,:)=[11500;98;0.2517];
SI_position(11,:)=[12000;93;0.709097];
SI_position(12,:)=[12500;91;0.8249];
SI_position(13,:)=[13000;92;0.7584];




%% 2. 数据重组:运行一次即可终身享用
% %将SI最后1s全部移到末端
% for k=1:length(SI_position)
%     DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(SI_position(k,1))];
%     DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
% end
% 
% for k=1:13
%     Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
%     xuhao_initial=[SI_position(k,2) SI_position(k,3)];
%     for  i_file=90:-1:1
%        xuhao_initial(1)=xuhao_initial(1)-1;
%        xuhao_initial(2)=1-xuhao_initial(2);
%        sst1=importdata(fullfile(DataInfo.location{k},['Compressor2Stall-',num2str(SI_position(k,1)),'-',num2str(xuhao_initial(1)),'.mat']));
%        sst2=importdata(fullfile(DataInfo.location{k},['Compressor2Stall-',num2str(SI_position(k,1)),'-',num2str(xuhao_initial(1)+1),'.mat']));
%        temp=sst1(round(xuhao_initial(2)*Process.sensorParameter.fs  ):end,:);
%        DATA=[temp;sst2(1:1*Process.sensorParameter.fs  -length(temp),:)];
%        save_directory = ['DATA/Compressor2Stall-',num2str(SI_position(k,1))];  %频谱图存储文件夹
%        save([save_directory,'/Transform-',num2str(SI_position(k,1)),'-',num2str(i_file),'.mat'],'DATA');
%     end
% end


%% 3. 计算阀门开度
%通过观察，阀门13000rpm在SI的阀门开度为29；7000rpm在SI的阀门开度为25;该值需要根据实验确定。
alpha_13000_t0=95;
alpha_13000_t90=29;
alpha_7000_t90=25;
alpha_13000_t=@(t) alpha_13000_t0-t.*(alpha_13000_t0-alpha_13000_t90)/90;
alpha_segm=(alpha_13000_t90-alpha_7000_t90)/12;
alpha_RS_t90=@(RS) alpha_13000_t90-alpha_segm*(13000-RS)./500;
alpha_RS_t0=@(RS) alpha_13000_t0-alpha_segm*(13000-RS)./500;





%% 4. 计算I_RI


%% 5. 换算I_RI
