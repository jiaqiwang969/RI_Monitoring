function [RotorSpeed]=waveletPlot_universal(signal,fs,RotorSpeed)    
sst1=signal;
variance = std(sst1)^2;
sst = (sst1 - mean(sst1))/sqrt(variance) ;
n = length(sst);
dt = 1/fs ;
round=[0:length(sst)-1]*dt*RotorSpeed/60;  % construct round array
xlim = [0,floor(round(end))];              % plotting range
pad = 1;                                   % pad the round series with zeroes (recommended)
dj = 0.125;                                % this will do 4 sub-octaves per octave
s0 = 15*dt/(fs/fs);                        % this says start at a scale of 6 months
j1 = 9/dj;                                 % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
mother = 'morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum
fk=1./period/(RotorSpeed/60);

for k_part=1:floor(round(end)/100)
    part{k_part}=find(round<100*k_part&round>5*(k_part-1));
end
for k_part=1:floor(round(end)/100)
    global_ws{k_part} =(sum(power(:,part{k_part})')./length(part{k_part}));   % round-average over all rounds
end


% fig=subplot('position',[0.77 0.11 0.2 0.53])
for k_part=[1]
plot((fk),global_ws{k_part},'DisplayName',[num2str((k_part-1)*5),'-',num2str((k_part)*5),'round'],'LineWidth',1,'Color',[1 0 0]); hold on
end

ylabel('Norm. Realtive Power (¡ª)')
title('c) Global Wavelet Spectrum')
Xticks = (fix((min(fk))):fix((max(fk))));
set(gca,'XLim',([min(fk),max(fk)]), ...
	'XDir','normal', ...
	'XTick',(Xticks(:)), ...
	'XTickLabel',Xticks)
% set(gca,'YLim',[0,1.25])
grid on;set(gca,'gridlinestyle','--');
hold on
% colorbar('peer',gca);
% legend1 = legend(gca,'off');
% set(legend1,...
%     'Position',[0.76926863456744 0.74909420289855 0.13115330636083 0.0612922705314016]);

