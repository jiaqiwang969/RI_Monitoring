clear all;
% 在一段正旋信号中间加入一段长度不同，且频率不同的正旋信号
t=0:0.001:4;
t_short=0:0.001:0.5;
s1=sin(2*pi*t);
s2=sin(20*pi*t);
s2_short=sin(20*pi*t_short);
s3=sin(2*pi*t);
s=[s1,s2,s3];    % 在一段正旋信号中间加入一段长度较长，且频率不同的正旋信号
s_short=[s1,s2_short,s3];    % 在一段正旋信号中间加入一段长度较短，且频率不同的正旋信号
figure(1);
subplot(211);
plot(s);
title('加入一段长度较长的不同频率正旋信号的原始信号');
s_fft=abs(fft(s));
s_ff=(0:round(length(s)/10))/(length(s)*(t(2)-t(1)));
subplot(212);
plot(s_ff,s_fft(1:round(length(s)/10)+1));
title('上面信号的傅立叶变换');
xlabel('傅立叶变换可以检测到两个频率');
figure(2);
subplot(211);
plot(s_short);
title('加入一段长度较短的不同频率正旋信号的原始信号');
s_short_fft=abs(fft(s_short));
s_short_ff=(0:round(length(s_short)/10))/(length(s_short)*(t(2)-t(1)));
subplot(212);
plot(s_short_ff,s_short_fft(1:round(length(s_short)/10)+1));
title('上面信号的傅立叶变换');
xlabel('傅立叶变换只可以检测到一个频率');
figure(3);
subplot(421);plot(s_short);
title('较短的原始信号');
ylabel('s-short');
[c,l]=wavedec(s_short,6,'db3');  %采用db3小波并对较短的原始信号进行六层分解
apcmp=wrcoef('a',c,l,'db3',6);
subplot(422);plot(apcmp);
ylabel('ca6');
for i=1:6
    decmp=wrcoef('d',c,l,'db3',7-i);
    subplot(4,2,2+i);
    plot(decmp);
    ylabel(['d',num2str(7-i)]);
end