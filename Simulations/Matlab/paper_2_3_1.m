%% 2.3.1 波片延迟量在波长、波数域的形貌及对应的通道分布
clc
clear
close all
addpath(genpath(pwd))
%% 基础参数设置
% 光谱仪分辨率
wvl_min = 450; % 最短波长，单位：nm
wvl_max = 800; % 最长波长，单位：nm
R_lambda = 0.5; % 波长分辨率，单位：nm
RP_lambda = ceil(( wvl_max - wvl_min )/R_lambda); % 波长分辨能力

% 计算系统波数域分辨率
wvn_min = 1/wvl_max; % 最小波数，单位：nm^-1
wvn_max = 1/wvl_min; % 最大波数，单位：nm^-1
wvn_range = wvn_max - wvn_min; % 有效波数范围，单位：nm^-1
R_sigma_max = R_lambda/(wvl_min^2); % 波数分辨率，最大值，单位：nm^-1
R_sigma_min = R_lambda/(wvl_max^2); % 波数分辨率，最小值，单位：nm^-1
R_sigma_avg = (R_sigma_max*R_sigma_min)^(1/2); % 波数分辨率，几何均值，单位：nm^-1
RP_sigma = ceil( ( wvn_max - wvn_min )/R_sigma_avg ); % 波数分辨能力

% 设置采样波数
wvn = wvn_min:( wvn_max - wvn_min )/(RP_sigma-1):wvn_max ;  % 系统采样波数，单位：nm^-1
wvl = flip(1./wvn);% 系统采样波长，单位：nm

% OPD域的可解析范围
OPD_max = (1/R_sigma_avg)/2; % OPD域有效的OPD范围，单位：nm
OPD = -OPD_max:2*OPD_max/(RP_sigma-1):OPD_max;
num = length(wvl);
win = blackman(num)';

%% 读取实测延迟量 
dataTable = readtable("Quartz_1.5mm_degree.csv","ReadVariableNames",1,"VariableNamingRule","preserve");
raw_wvl = dataTable.lambda;
% 延迟量单位：度
delta_t24 = spline(raw_wvl,dataTable.t24,wvl);

figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 创建第一个坐标系（左下）
ax1 = axes('Position', [0.01 0.01 0.65*1.5 0.815*1.2]);
plot(ax1, wvl,delta_t24,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-');    % 左侧曲线
xlim(ax1, [430, 830]); 
ylim([5000,14000])
box(ax1, 'off');           % 关闭边框
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[]);

% 创建第二个坐标系（右上）
ax2 = axes('Position', [0.01 0.01 0.65*1.5 0.815*1.2]);
plot(ax2, wvn*1000,flip(delta_t24),'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','--');    % 右侧曲线
xlim(ax2, [1.2, 2.3]);        % 设置右侧 x 轴范围
box(ax2, 'off');           % 关闭边框

% 设置第二个坐标系透明背景
set(ax2, 'Color', 'none', 'XAxisLocation', 'top', 'YAxisLocation', 'right');

ylim([5000,14000])
% 添加图例
leg = legend([ax1.Children, ax2.Children],'波长域','波数域');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2;
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[]);
% print('fig2_3_a','-dpng','-r600')

cos_d = cosd(delta_t24);
figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001,abs(fftshift(fft(cos_d.*win))),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75, 75])
ylim([-10,160])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig2_3_b','-dpng','-r600')


