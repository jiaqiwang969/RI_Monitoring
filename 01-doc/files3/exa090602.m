% -------------------------------------------------------------------
% exa060602.m, for fig 9.6.2
% ���� morlet С����ʱ���μ�Ƶ�� 
%--------------------------------------------------------------------
clear;

lb=-4;ub=4;n=1000;
[wavelet,x]=morlet(lb,ub,n);
subplot(2,2,1)
plot(x,wavelet);grid on
title(' Morlet wavelet: Psi')

a=[1];
fs=1/(x(2)-x(1));
[h,w]=freqz(wavelet,a,2048,'whole',fs);
hr=abs(h);
subplot(2,2,2)
plot(w(1:128)*2*pi,hr(1:128));grid on
title('The FT of Psi')


