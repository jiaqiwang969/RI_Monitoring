
function [Tdata_resample]=dataReconstruct(Data,Rotor_Speed,Pulse,resamplePoint,fs)

H_if_RI=1:(length(Rotor_Speed)-1);
solutime=1;
for H=H_if_RI  %ÿתһȦˢ��һ������һ��ʱ�䲽��%����RI_marker����
%% ����Ͱ������λ�ù�ϵ
distance_key_rotor=6;%����λ�ú����ҶƬ֮��ľ���
R_rotor=184.4002;%mm/��Ҷ�뾶
x=[2;4.41;7.68;11.14;14.39;17.41;20.73;23.77;26.63;34];
y_2=[20.29;24.92;29.60;35.37;40.80;44.83;49.71;55.40;59.70;72.36];
distance_point_B1_R1=3.8*R_rotor*2*pi/29;
distance_point_B1_C1=5.3/360*R_rotor*2*pi;%17.0575mm;��C1�Ƶ�70mm����ֵΪ��70-17.0575��mm
distance_pointR1_R1=R_rotor*2*pi/29;
interval=Rotor_Speed(H)/60*R_rotor/fs*2*pi;
round_point=distance_pointR1_R1*29/interval;%һ�ܵĵ���
Len=floor(round_point)-2;%=130;%������ʾ�ڵ����;ȡһ��Χ��һ����ά��Ȧ
point_Pulse_rotor=round(distance_key_rotor/interval); %%%12000-��19*1.6092����13000-��12*1.7433����10000-��5*1.0740����11000-��-2*1.3598����125000-��16*1.6763��
point_rotor_rotor=round_point/29; %����ҶƬ֮��ľ������������һ��Ҷ���ĵ��� %144000/BPF
tran_C1=(70-distance_point_B1_C1)/interval;%C1ʵ�ʵ�λ�������λ��23.2235����λ��C1ǰ��23����λ,��24
point_B1_R1_interval=round(distance_point_B1_R1./interval);  %B1��R1��52.5761����λ��R1-R8ǰ��53����λ����54
point_B1_C1_interval=round((y_2(end)-y_2(1)-distance_point_B1_C1)./interval);  %%17.0575mm;��C1�Ƶ�70mm����ֵΪ��70-17.0575��mm
%% �����������
y_1=(1:Len+1)*interval;
X = repmat(x,size(y_1))';
yy_1 = repmat(y_1,size(x))';
yy_2=repmat(y_2,size(y_1))';
Y=yy_1+yy_2;
% 2ά���3ά��ת��
initial_y=Y(1,1);
Theta=(Y-initial_y)./R_rotor;
XX=R_rotor*cos(Theta);
YY=R_rotor*sin(Theta);
ZZ=X;

% �������ݣ�������ֵ
    rpm{solutime}(:,1) = Data(Pulse(H):Pulse(H)+Len,4);%������B1
    rpm{solutime}(:,2) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),8);
    rpm{solutime}(:,3) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),9);
    rpm{solutime}(:,4) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),10);
    rpm{solutime}(:,5) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),11);
    rpm{solutime}(:,6) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),12);
    rpm{solutime}(:,7) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),13);
    rpm{solutime}(:,8) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),14);
    rpm{solutime}(:,9) = Data((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),15);
    rpm{solutime}(:,10) = Data((point_B1_C1_interval+Pulse(H)):(point_B1_C1_interval+Pulse(H)+Len),16);
   % һȦ��29��ҶƬ��������л��֣���Ƕ�Ҷλ�ã�29�ݣ�:���
      for k=1:29  
        Pulse_solid(k)=round(1+point_rotor_rotor*(k-1));%�����Ϳ����ڱ�ת�ٲ���ҶƬ��
      end
     rpm{solutime}(Pulse_solid+point_Pulse_rotor,:) =1.1*max(rpm{solutime}(:,end));  %�ڴ�Լpoint_rotor_rotor������Ѱ��ȷ��ҶƬλ�ã��ݶ��÷�����Ϊ����ָ��

 %%
     tdata.surfaces(solutime).zonename='mysurface zone';
     tdata.surfaces(solutime).x=XX;    %size 3x3 
     tdata.surfaces(solutime).y=ZZ;    %size 3x3
     tdata.surfaces(solutime).z=YY;    %size 3x3
     tdata.surfaces(solutime).v(1,:,:)=rpm{solutime};%���ݼ����ź��ж�
     tdata.surfaces(solutime).solutiontime=solutime;
     Tdata_resample.surfaces(solutime).v=resample(rpm{solutime},resamplePoint,length(rpm{solutime}));%���ݼ����ź��ж�
     %Tdata_nonresample.surfaces(solutime).v=rpm{solutime};
     solutime=solutime+1;
end

end