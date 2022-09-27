
function [Tdata_resample]=dataReconstruct(Data,Rotor_Speed,Pulse,resamplePoint,fs)

H_if_RI=1:(length(Rotor_Speed)-1);
solutime=1;
for H=H_if_RI  %每转一圈刷新一次网格，一个时间步！%按照RI_marker排序
%% 计算和安排相对位置关系
distance_key_rotor=6;%键向位置和最近叶片之间的距离
R_rotor=184.4002;%mm/动叶半径
x=[2;4.41;7.68;11.14;14.39;17.41;20.73;23.77;26.63;34];
y_2=[20.29;24.92;29.60;35.37;40.80;44.83;49.71;55.40;59.70;72.36];
distance_point_B1_R1=3.8*R_rotor*2*pi/29;
distance_point_B1_C1=5.3/360*R_rotor*2*pi;%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm
distance_pointR1_R1=R_rotor*2*pi/29;
interval=Rotor_Speed(H)/60*R_rotor/fs*2*pi;
round_point=distance_pointR1_R1*29/interval;%一周的点数
Len=floor(round_point)-2;%=130;%调节显示节点个数;取一周围成一个三维的圈
point_Pulse_rotor=round(distance_key_rotor/interval); %%%12000-（19*1.6092）；13000-（12*1.7433）；10000-（5*1.0740）；11000-（-2*1.3598）；125000-（16*1.6763）
point_rotor_rotor=round_point/29; %两个叶片之间的距离来推算采样一个叶道的点数 %144000/BPF
tran_C1=(70-distance_point_B1_C1)/interval;%C1实际点位和理想点位差23.2235个点位，C1前移23个点位,第24
point_B1_R1_interval=round(distance_point_B1_R1./interval);  %B1与R1差52.5761个点位，R1-R8前移53个点位，第54
point_B1_C1_interval=round((y_2(end)-y_2(1)-distance_point_B1_C1)./interval);  %%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm
%% 生成网格矩阵
y_1=(1:Len+1)*interval;
X = repmat(x,size(y_1))';
yy_1 = repmat(y_1,size(x))';
yy_2=repmat(y_2,size(y_1))';
Y=yy_1+yy_2;
% 2维面和3维面转换
initial_y=Y(1,1);
Theta=(Y-initial_y)./R_rotor;
XX=R_rotor*cos(Theta);
YY=R_rotor*sin(Theta);
ZZ=X;

% 导入数据，对网格赋值
    rpm{solutime}(:,1) = Data(Pulse(H):Pulse(H)+Len,4);%传感器B1
    rpm{solutime}(:,2) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),8);
    rpm{solutime}(:,3) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),9);
    rpm{solutime}(:,4) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),10);
    rpm{solutime}(:,5) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),11);
    rpm{solutime}(:,6) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),12);
    rpm{solutime}(:,7) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),13);
    rpm{solutime}(:,8) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),14);
    rpm{solutime}(:,9) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),15);
    rpm{solutime}(:,10) = Data((point_B1_C1_interval+Pulse(H)):(point_B1_C1_interval+Pulse(H)+Len),16);
   % 一圈有29个叶片，对其进行划分，标记动叶位置（29份）:标记
      for k=1:29  
        Pulse_solid(k)=round(1+point_rotor_rotor*(k-1));%这样就可以在变转速插入叶片了
      end
     rpm{solutime}(Pulse_solid+point_Pulse_rotor,:) =1.1*max(rpm{solutime}(:,end));  %在大约point_rotor_rotor个点搜寻正确的叶片位置，暂定用方差作为衡量指标

 %%
     tdata.surfaces(solutime).zonename='mysurface zone';
     tdata.surfaces(solutime).x=XX;    %size 3x3 
     tdata.surfaces(solutime).y=ZZ;    %size 3x3
     tdata.surfaces(solutime).z=YY;    %size 3x3
     tdata.surfaces(solutime).v(1,:,:)=rpm{solutime};%根据键向信号判读
     tdata.surfaces(solutime).solutiontime=solutime;
     Tdata_resample.surfaces(solutime).v=resample(rpm{solutime},resamplePoint,length(rpm{solutime}));%根据键向信号判读
     %Tdata_nonresample.surfaces(solutime).v=rpm{solutime};
     solutime=solutime+1;
end

end