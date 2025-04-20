%% 4.4.1.3 研究随机噪声对装调方法的影响
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

% ch0 = [344, 358];
ch2 = [373, 387];
ch3 = [388, 401];
ch4 = [402, 416];

seed = 5;
rng(seed);

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

%% 方位角变化量与光谱强度变化量的关系
% 重置系统方位角--调整多级波片1
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 0;
Iout = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
ch = abs(fftshift(fft(Iout.*win)));
ch3_ideal_value = max(ch(ch3(1):ch3(2)));
plot(abs(fftshift(fft(Iout.*win))));
for j = 1:10
    % 由于是随机噪声，重复十次取平均值
    for i=1:9
        SNR = 10+5*(i-1);
        Iout_2 = awgn(Iout, SNR, 'measured');
        ch = abs(fftshift(fft(Iout_2.*win)));
        ch3_delta_value(i,j) = (max(ch(ch3(1):ch3(2)))-ch3_ideal_value)/ch3_ideal_value;
    end
end
figure(1)
plot(10:5:50, abs(mean(ch3_delta_value,2)))

% 重置系统方位角--调整检偏器
configuration.P_azimuth = 0;
configuration.r1_azimuth = 45;
configuration.r2_azimuth = 0;
configuration.A_azimuth = 45;
Iout = instruModel(phi1,phi2,Iin,M_air,0,num,configuration);
ch = abs(fftshift(fft(Iout.*win)));
ch2_ideal_value = max(ch(ch2(1):ch2(2)));
ch3_ideal_value = max(ch(ch3(1):ch3(2)));
ch4_ideal_value = max(ch(ch4(1):ch4(2)));
for j = 1:10
    for i=1:9
        SNR = 10+5*(i-1);
        Iout_2 = awgn(Iout, SNR, 'measured');
        ch = abs(fftshift(fft(Iout_2.*win)));
        ch2_delta_value(i,j) = (max(ch(ch2(1):ch2(2)))-ch2_ideal_value)/ch2_ideal_value;
        ch4_delta_value(i,j) = (max(ch(ch4(1):ch4(2)))-ch4_ideal_value)/ch4_ideal_value;
    end
end
figure(2)
plot(10:5:50, abs(mean(ch2_delta_value,2)), 10:5:50, abs(mean(ch4_delta_value,2)))
legend('ch2','ch4')

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(10:5:50, abs(mean(ch3_delta_value,2))*100,'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([5,55])
ylim([-0.1,1.35])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_10_a','-dpng','-r600')

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(10:5:50, abs(mean(ch2_delta_value,2))*100,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-');
hold on
plot(10:5:50, abs(mean(ch4_delta_value,2))*100,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
xlim([5,55])
ylim([-0.1,2.2])
leg = legend("±2频", "±4频");
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 1; leg.Location = "best";
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","宋体", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_4_10_b','-dpng','-r600')
