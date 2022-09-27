clc
clear
% Detection of Transients
rng default;
dt = 0.001;
t = 0:dt:1-dt;
N = length(t);
fs = 1/dt;
addNoise = 0.025*randn(size(t));
x = cos(2*pi*150*t).*(t>=0.1 & t<0.3)+sin(2*pi*200*t).*(t>0.7);
x = x+addNoise;
x([222 800]) = x([222 800 ])+[-2 2];
figure;
plot(t.*1000,x);
xlabel('Milliseconds'); ylabel('Amplitude');

figure;
[cfs,f_cwt] = cwt(x,fs,'amor');
contour(t*1000,f_cwt,abs(cfs))
grid on
c = colorbar;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt');
c.Label.String = 'Magnitude';

figure
surf(t.*1000,f_cwt,abs(cfs),'edgecolor','none'); view(0,90);
axis tight;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt');

% figure
% imagesc(t.*1000,f_cwt,abs(cfs));
% axis xy;
% xlabel('Milliseconds')
% ylabel('Frequency')
% title('cwt');

figure;
f = linspace(1, fs/2, N/2);
Spec_cwt = cwt_cmor(x,1,2,f,fs);
contour(t*1000,f,abs(Spec_cwt))
grid on
c = colorbar;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt\_cmor');
c.Label.String = 'Magnitude';

figure
surf(t.*1000,f,abs(Spec_cwt),'edgecolor','none'); view(0,90);
axis tight;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt\_cmor');

figure;
f = linspace(1, fs/2, N/2);
Spec_cwt = cwt_cmor_2(x,1,2,f,fs);
contour(t*1000,f,abs(Spec_cwt))
grid on
c = colorbar;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt\_cmor2');
c.Label.String = 'Magnitude';

figure
surf(t.*1000,f,abs(Spec_cwt),'edgecolor','none'); view(0,90);
axis tight;
xlabel('Milliseconds')
ylabel('Frequency')
title('cwt\_cmor2');


% figure
% imagesc(t.*1000,f,abs(Spec_cwt).^2); 
% axis xy;
% xlabel('Milliseconds')
% ylabel('Frequency')
% title('cwt\_cmor');

figure
[S,F,T] = spectrogram(x,50,48,128,1000);
surf(T.*1000,F,abs(S),'edgecolor','none'); view(0,90);
axis tight;
xlabel('Milliseconds')
ylabel('Frequency')
title('STFT')

figure
contour(T.*1000,F,abs(S))
axis tight;
xlabel('Time'), ylabel('Hz');
title('STFT')

