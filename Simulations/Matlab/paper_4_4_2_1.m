%% 4.4.2.1 演示延迟量的校准过程
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

ch0 = [348, 354];
ch2 = [373, 387];
ch3 = [388, 401];
ch4 = [402, 416];

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

%% 进行波片延迟量求解
% 重置系统方位角--调整多级波片1
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
Iout_1 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
f_1 = abs(fftshift(fft(Iout_1.*win)));
% figure(1)
% plot(wvn, Iout_1)
% figure(2)
% plot(f_1)

% 重置系统方位角--调整检偏器
configuration.A_azimuth = 45;
Iout_2 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
f_2 = abs(fftshift(fft(Iout_2.*win)));
% figure(3)
% plot(wvn, Iout_2)
% figure(4)
% plot(OPD, f_2)

% 进行延迟量求解
[e_delta1, e_delta1minus2, e_delta1plus2, e_delta2] = channel_calibration_auto(Iout_1, Iout_2, ch0, ch2, ch3, ch4);
% 校准延迟量
d1_cali = unwrap(angle(e_delta1)); % rad
d1_cali = -d1_cali*180/pi;
d2_cali = unwrap(angle(e_delta2)); % rad
d2_cali = -d2_cali*180/pi;
cos_d1 = cosd(d1_cali); sin_d1 = sind(d1_cali);
cos_d2 = cosd(d2_cali); sin_d2 = sind(d2_cali);

% 对比校准量和真实值
% figure(5)
% plot(wvn, cos_phi_1, '--', wvn, cos_d1, '-');
% xlim([1.2, 2.3]*1e-3)
% ylim([-1.2, 2])
% legend('cosphi_1', 'cosphi_1cali')
% figure(6)
% plot(wvn, cos_phi_2, '--', wvn, cos_d2, '-');
% xlim([1.2, 2.3]*1e-3)
% ylim([-1.2, 2])
% legend('cosphi_2', 'cosphi_2cali')

% 可视化RMSE随波段范围的变化
for i = 1:70
    match_range = 1+i:num-i;
    RMSE_1(i) = RMSE_single(cos_d1, cos_phi_1, match_range);
end

for i = 1:70
    match_range = 1+i:num-i;
    RMSE_2(i) = RMSE_single(cos_d2, cos_phi_2, match_range);
end
% figure(7)
% plot(1:70, RMSE_1, 1:70, RMSE_2)
% legend('phi_1','phi_2')

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,4.5/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000, Iout_1,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([1.2,2.3])
ylim([-0.1,0.6])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_a_1','-dpng','-r600')

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,4.5/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001, f_1,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75,75])
ylim([-5,90])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_a_2','-dpng','-r600')

figure(3)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000, cos_phi_1,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','--');
hold on
plot(wvn*1000, cos_d1,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([1.2,2.3])
ylim([-1.2,1.5])
leg = legend("真实值", "校准值");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_b','-dpng','-r600')

figure(4)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,4.5/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000, Iout_2,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([1.2,2.3])
ylim([-0.1,0.6])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_c_1','-dpng','-r600')

figure(5)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,4.5/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(OPD*0.001, f_2,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-75,75])
ylim([-5,90])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_c_2','-dpng','-r600')

figure(6)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn*1000, cos_phi_2,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','--');
hold on
plot(wvn*1000, cos_d2,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([1.2,2.3])
ylim([-1.2,1.5])
leg = legend("真实值", "校准值");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_11_d','-dpng','-r600')

figure(7)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(1:70, RMSE_1,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','--');
hold on
plot(1:70, RMSE_2,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([0,71])
ylim([0,0.25])
leg = legend("R_1 校准误差", "R_2 校准误差");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_12','-dpng','-r600')