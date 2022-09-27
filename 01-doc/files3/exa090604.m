% -------------------------------------------------------------------
% exa060604.m, for fig 9.6.4
% 产生 gauswavf 小波的时域波形及频谱 
%--------------------------------------------------------------------
clear;
figure
[f1,x1]=gauswavf(-6,6,128,2);
[f2,x2]=gauswavf(-6,6,128,6);
subplot(221)
%stem(f)
plot(x1,f1);grid on
hold on
plot(x2,f2)
title('Gaussian wavelet: Psi')
a=[1];
   fs=1;
[h1,w1]=freqz(f1,a,256,'whole',fs);
hr1=abs(h1);
[h2,w2]=freqz(f2,a,256,'whole',fs);
hr2=abs(h2);
subplot(222)
plot(w1(1:128),hr1(1:128));grid on
hold on
plot(w2(1:128),hr2(1:128))
title(' The FT of Psi')


