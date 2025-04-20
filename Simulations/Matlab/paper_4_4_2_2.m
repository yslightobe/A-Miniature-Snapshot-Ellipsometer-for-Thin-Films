%% 4.4.2.2 研究校准过程对随机噪声的鲁棒性
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

% 重置系统方位角--调整检偏器
configuration.A_azimuth = 45;
Iout_2 = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
f_2 = abs(fftshift(fft(Iout_2.*win)));

% 可视化RMSE随信噪比的变化
match_range = 50:num-50;
for j = 1:10
    for i=1:9
        SNR = 10+5*(i-1);
        Iout_1_error = awgn(Iout_1, SNR, 'measured');
        Iout_2_error = awgn(Iout_2, SNR, 'measured');
        % 进行延迟量求解
        [e_delta1, e_delta1minus2, e_delta1plus2, e_delta2] = channel_calibration_auto(Iout_1_error, Iout_2_error, ch0, ch2, ch3, ch4);
        % 校准延迟量
        d1_cali = unwrap(angle(e_delta1)); % rad
        d1_cali = -d1_cali*180/pi;
        d2_cali = unwrap(angle(e_delta2)); % rad
        d2_cali = -d2_cali*180/pi;
        cos_d1 = cosd(d1_cali);
        cos_d2 = cosd(d2_cali);
        % 计算RMSE
        RMSE_1(i, j) = RMSE_single(cos_d1, cos_phi_1, match_range);
        RMSE_2(i, j) = RMSE_single(cos_d2, cos_phi_2, match_range);
    end
end
figure(1)
plot(10:5:50, abs(mean(RMSE_1,2)), 10:5:50, abs(mean(RMSE_2,2)))
legend('phi_1', 'phi_2')

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(10:5:50, abs(mean(RMSE_1,2)),'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','--');
hold on
plot(10:5:50, abs(mean(RMSE_2,2)),'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([5,55])
ylim([0,0.25])
leg = legend("R_1 校准误差", "R_2 校准误差");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_13','-dpng','-r600')