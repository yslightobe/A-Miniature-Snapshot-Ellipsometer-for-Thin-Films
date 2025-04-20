function createfigure(wavenumber,N,N_RC2,C,C_RC2,S,S_RC2,match_range,...
    legOrnot, absolute_error_range, xlabelOrnot, ylim_1)

figure
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16,9]);
subplot(1,2,1)
% 使用 plot 的矩阵输入创建多行
plot(wavenumber,N,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(wavenumber,C,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-'); hold on
plot(wavenumber,S,'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-'); hold on
plot(wavenumber,N_RC2,'LineWidth',0.8,'Color',[0 0 0],'LineStyle','--'); hold on
plot(wavenumber,C_RC2,'LineWidth',0.8,'Color',[0 0 0],'LineStyle','--'); hold on
plot(wavenumber,S_RC2,'LineWidth',0.8,'Color',[0 0 0],'LineStyle','--');

if legOrnot
    leg = legend('{\itN}','{\itC}','{\itS}','truth');
    leg.ItemTokenSize = [6,5]; leg.Box = 'off'; leg.NumColumns = 4;
end

xlim([1.2 2.3]);
% ylim([-1.2 1.2]);
% ylim([-1 1.5]);
ylim(ylim_1);

if xlabelOrnot
    set(gca,'unit','normalized','InnerPosition',[0.16,0.6,0.28,0.3],'FontName',...
        'Times New Roman','FontSize',12,'LineWidth',1,'XTick',...
        [1.4 1.6 1.8 2 2.2 2.4],'XTickLabel',...
        {'1.4','1.6','1.8','2','2.2','2.4'},'Box','on','XTickLabelRotation',0);
else
    set(gca,'unit','normalized','InnerPosition',[0.16,0.6,0.28,0.3],'FontName',...
        'Times New Roman','FontSize',12,'LineWidth',1,'XTick',...
        [1.4 1.6 1.8 2 2.2 2.4],'XTickLabel',...
        {'','','','','',''},'Box','on','XTickLabelRotation',0);
end

subplot(1,2,2)
boxplot([abs(N(match_range)-N_RC2(match_range)) abs(C(match_range)-C_RC2(match_range)) ...
    abs(S(match_range)-S_RC2(match_range))],...
     'Symbol','', 'ColorGroup',[0.741,0,0.447])

ylim(absolute_error_range);

set(gca,'unit','normalized','InnerPosition',[0.53,0.6,0.28,0.3],...
    'FontName','Times New Roman','FontSize',12,'LineWidth',1,...
    'TickLabelInterpreter','none','XTick',[],'XTickLabel',{'','',''});