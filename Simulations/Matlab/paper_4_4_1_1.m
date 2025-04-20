%% 4.4.1.1 研究装调过程中通道分布的变化
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
Iin = ones(1, num);
win = blackman(num)';

configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 45;
M_air = repmat(eye(4), [1,1,num]);

%% 读取数据
% 多级波片延迟量
dataTable = readtable("Quartz_1.5mm_degree.csv","ReadVariableNames",1,"VariableNamingRule","preserve");
raw_wvl = dataTable.lambda;
phi2 = spline(raw_wvl,dataTable.t24,wvl);
phi1 = 3*phi2;

cos_phi_2 = cosd(phi2); sin_phi_2 = sind(phi2);
cos_phi_1 = cosd(phi1); sin_phi_1 = sind(phi1);
cos_phi1_s2 = cos_phi_1.*sin_phi_2;
sin_phi1_s2 = sin_phi_1.*sin_phi_2;

%% 安装两个多级波片时通道分布的变化
configuration.P_azimuth = 0;
configuration.r1_azimuth = 30;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
Iout = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
configuration.P_azimuth = 0;
configuration.r1_azimuth = 0;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
Iout_0 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
configuration.P_azimuth = 0;
configuration.r1_azimuth = 0;
configuration.r2_azimuth = 30;
configuration.A_azimuth = 0;
Iout_1 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
% % 0-0-0-0
% figure(1)
% plot(abs(fftshift(fft(Iout_0.*win))))
% % 0-30-0-0
% figure(2)
% plot(abs(fftshift(fft(Iout.*win))))
% % 0-0-30-0
% figure(3)
% plot(abs(fftshift(fft(Iout_1.*win))))

%% 转动方位角时，对应频道强度变化
% 0-0-0-0 → 0-45-0-0
configuration.P_azimuth = 0;
configuration.r1_azimuth = 0;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
xx = 50;
for i=1:xx
    configuration.r1_azimuth = 0+45*(i-1)/(xx-1);
    Iout_r1 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
    f_r1 = abs(fftshift(fft(Iout_r1.*win)));
    ch0_r1(i) = max(f_r1(348:354));
    ch1_r1(i) = max(f_r1(360:372));
    ch2_r1(i) = max(f_r1(375:385));
    ch3_r1(i) = max(f_r1(388:401));
    ch4_r1(i) = max(f_r1(402:415));
end
% 0-45-0-0 → 0-45-0-45
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
xx = 50;
for i=1:xx
    configuration.A_azimuth = 0+45*(i-1)/(xx-1);
    Iout_a = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
    f_a = abs(fftshift(fft(Iout_a.*win)));
    ch0_a(i) = max(f_a(348:354));
    ch1_a(i) = max(f_a(360:372));
    ch2_a(i) = max(f_a(375:385));
    ch3_a(i) = max(f_a(388:401));
    ch4_a(i) = max(f_a(402:415));
end
figure(4)
plot(1:xx, ch3_r1)
figure(5)
plot(1:xx, ch2_a, 1:xx, ch3_a, 1:xx, ch4_a)
legend("2频","3频","4频")

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,8/2.3,6/2.3]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001, abs(fftshift(fft(Iout_0.*win))),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75,75])
ylim([-10, 170])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_6_b2','-dpng','-r600')

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,8/2.3,6/2.3]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001, abs(fftshift(fft(Iout.*win))),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75,75])
ylim([-10, 170])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_6_b1','-dpng','-r600')

figure(3)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,8/2.3,6/2.3]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001, abs(fftshift(fft(Iout_1.*win))),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75,75])
ylim([-10, 170])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_6_d1','-dpng','-r600')

figure(4)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(1:xx, ch0_r1,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-');
hold on
plot(1:xx, ch1_r1,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
plot(1:xx, ch2_r1,'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-');
plot(1:xx, ch3_r1,'LineWidth',1.2,'Color',[25 153 178]/255,'LineStyle','-');
plot(1:xx, ch4_r1,'LineWidth',1.2,'Color',[0.580, 0, 0.827],'LineStyle','-');
xlim([0,51])
ylim([-10,200])
% leg = legend('±0','±1','±2','±3','±4');
% leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 5;
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_7_a','-dpng','-r600')

figure(5)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(1:xx, ch0_a,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-');
hold on
plot(1:xx, ch1_a,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
plot(1:xx, ch2_a,'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-');
plot(1:xx, ch3_a,'LineWidth',1.2,'Color',[25 153 178]/255,'LineStyle','-');
plot(1:xx, ch4_a,'LineWidth',1.2,'Color',[0.580, 0, 0.827],'LineStyle','-');
xlim([0,51])
ylim([-10,200])
leg = legend('0','±1','±2','±3','±4');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 5;
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_7_b','-dpng','-r600')