%% 3.2.1 不同厚度SiO2的参数和石英波片延迟量变化
% 无方位角误差和温度波动，二氧化硅薄膜厚度为0.03um
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

% 生成不同厚度薄膜的参数
[~, N_10, ~, ~] = rcwa_film(wvl/1000,0.01,65); % 微米(行向量)，微米，角度
[~, N_100, ~, ~] = rcwa_film(wvl/1000,0.1,65);
[~, N_1000, ~, ~] = rcwa_film(wvl/1000,1,65);

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000,N_10,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(wvn*1000,N_100,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
plot(wvn*1000,N_1000,'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-');
ylim([-1,2])
leg = legend('10nm','100nm','1000nm');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 3;

set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig3_1_a','-dpng','-r600')

%% 读取实测延迟量 
dataTable = readtable("Quartz_1.5mm_degree.csv","ReadVariableNames",1,"VariableNamingRule","preserve");
raw_wvl = dataTable.lambda;
% 延迟量单位：度
delta_t20 = spline(raw_wvl,dataTable.t20,wvl);
delta_t24 = spline(raw_wvl,dataTable.t24,wvl);
delta_t28 = spline(raw_wvl,dataTable.t28,wvl);
delta_t32 = spline(raw_wvl,dataTable.t32,wvl);
delta_t35 = spline(raw_wvl,dataTable.t35,wvl);

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000,delta_t20-delta_t24,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(wvn*1000,delta_t28-delta_t24,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
plot(wvn*1000,delta_t32-delta_t24,'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-');
plot(wvn*1000,delta_t35-delta_t24,'LineWidth',1.2,'Color',[25 153 178]/255,'LineStyle','-');
ylim([-20,15])
leg = legend('20℃','28℃','32℃','35℃');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 4;

set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig3_1_b','-dpng','-r600')