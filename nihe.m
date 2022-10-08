
jet_color2=colormap(jet(15));
A=[];
for k=3:18
    A(k,:)=UI_FFT_noAverage(k,DataInfo.condition(k,2)-89:DataInfo.condition(k,2)-89+78,2);
    plot(A(k,:),'Color',jet_color2(k,:),'LineWidth',2)
    hold on
end
%
% scale=1;
% plot([1:67]+14,1.4*smooth(UI_FFT(16,1:67,2)),'Color','k','LineWidth',2)
% plot([1:76]+4, 1.4*smooth(UI_FFT(17,1:76,2)),'Color','b','LineWidth',2)
% plot([1:75]+5, 1.4*smooth(UI_FFT(18,1:75,2)),'Color','y','LineWidth',2)

legend1 = legend('7000','7500','8000','8500','9000','9500','10000','10500','11000','11500','12000','12500','13000','14000','15000','16000')
set(legend1,'Location','northwest');
famen=90:-5:29;
xuhao=round(famen*DataInfo_a(1)+DataInfo_b(1))-10;
set(gca,'FontSize',14,'XGrid','on','XTick',xuhao,...
    'XTickLabel',mat2cell(famen.',ones(length(famen),1)));
grid on
    txt = {['转速','7000:500:13000','rpm'],...
        ['采样',num2str(N_caiyang*29),'点/圈'],['位置',DataInfo.sensorArray{2}]};
    annotation('textbox',...
        [0.3 0.77 0.17 0.13],...
        'String',txt,...
        'FontSize',13,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    set(gca,'FontSize',20);
% 创建 ylabel
% ylabel({'PI（）'});
ylabel('$I_{RI}=\sum_{i}P(i)\Delta f[pa^2]$','interpreter','latex','FontSize',20)
% 创建 xlabel
xlabel('磁阀开度(%)','FontSize',30);
% 设置其余坐标区属性
title('不同转速的I_R_I','FontSize',14)
set(gca,'FontSize',14);
saveas(h1,[save_directory,'/','FIG-所有转速-无平均-方差','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')

