%% 4.4.1.2 研究装调方法的灵敏度
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
%% 偏导数值求解及可视化
lamda_num = (700-450+1);
theta_R1 = 0:0.1:45;
theta_R1_2 = 90:0.1:135;
deriv_R1 = 0.5*(cos_phi_1(lamda_num)-1)*sind(4*theta_R1);
deriv_R1_2 = 0.5*(cos_phi_1(lamda_num)-1)*sind(4*theta_R1_2);

% figure(1)
% plot(theta_R1, abs(deriv_R1))
% figure(2)
% plot(theta_R1_2, abs(deriv_R1_2))

theta_A = 0:0.1:45;
deriv_A = 0.5*(-sin_phi_1(lamda_num)*sin_phi_2(lamda_num)*cosd(2*theta_A)-cos_phi_1(lamda_num)*sind(2*theta_A));
deriv_A_2 = 0.5*(sin_phi_1(lamda_num)*sin_phi_2(lamda_num)*cosd(2*theta_A)-cos_phi_1(lamda_num)*sind(2*theta_A));

% figure(3)
% plot(theta_A, abs(deriv_A), theta_A, abs(deriv_A_2))
% legend("取+号", "取-号")
%% 方位角变化量与光谱强度变化量的关系
% 重置系统方位角--调整多级波片1
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
Iout = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
for i = 1:10
    d_r1_azimuth = i*0.1;
    configuration.r1_azimuth = 45 - d_r1_azimuth;
    Iout_2 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
    delta_Iout_r1(i) = abs(Iout(lamda_num) - Iout_2(lamda_num));
end
figure(4)
plot(0.1:0.1:1, delta_Iout_r1)

% 重置系统方位角--调整检偏器A
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 45;
Iout = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
for i = 1:10
    d_a_azimuth = i*0.1;
    configuration.A_azimuth = 45 - d_a_azimuth;
    Iout_3 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
    delta_Iout_a(i) = abs(Iout(lamda_num) - Iout_3(lamda_num));
end
% figure(5)
% plot(0.1:0.1:1, delta_Iout_a)

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(theta_R1, abs(deriv_R1),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-1,46])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_8_a','-dpng','-r600')

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(0.1:0.1:1, delta_Iout_r1,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([0.,1.1])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_8_b','-dpng','-r600')

figure(3)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(theta_A, abs(deriv_A),'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(theta_A, abs(deriv_A_2),'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([-1,46])
ylim([0, 0.5])
leg = legend("取+号", "取-号");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 1; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_9_a','-dpng','-r600')

figure(4)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(0.1:0.1:1, delta_Iout_a,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([0.,1.1])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_9_b','-dpng','-r600')