clc
clear
[Lo_D,Hi_D]=wfilters('db7');
[Hi,Fq1] = freqz(Hi_D);
[Lo,Fq2] = freqz(Lo_D);
figure
plot(Fq1/2/pi,abs(Hi)/sqrt(2));
hold on
plot(Fq2/2/pi,abs(Lo)/sqrt(2));
set(gca,'Xlim',[0,0.5])