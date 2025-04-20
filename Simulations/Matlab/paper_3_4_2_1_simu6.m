%% 3.4.2.1 仿真实验#6 验证各算法对随机噪声的鲁棒性
% 方位角误差为1.5°、3.5°、1.3°、3.7°，无温度波动，二氧化硅薄膜厚度为0.8um
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

% 添加方位角误差
configuration.P_azimuth = 0+1.5;
configuration.r1_azimuth = 45+3.5;
configuration.r2_azimuth = 0+1.3;
configuration.A_azimuth = 45+3.7;

hyper_paras1 = [0.05, 0.02];
hyper_paras2 = [0.05, 0.02, 0.01, 1, 0]; % 最后两个参数控制是否进行方位角和温度补偿

%% 读取参考校准延迟量（24℃），单位：度
dataTable = readtable("Quartz_1.5mm_degree.csv","ReadVariableNames",1,"VariableNamingRule","preserve");
raw_wvl = dataTable.lambda;
phi2 = spline(raw_wvl,dataTable.t24,wvl);
phi1 = 3*phi2;
cos_phi_2 = cosd(phi2); sin_phi_2 = sind(phi2);
cos_phi_1 = cosd(phi1); sin_phi_1 = sind(phi1);
sin_phi1_s2 = sin_phi_1.*sin_phi_2;
sin_phi1_c2 = sin_phi_1.*cos_phi_2;

%% 读取测量时的延迟量（≠24℃），单位：度 此时没有温度波动，取24℃
phi2_realtime = spline(raw_wvl,dataTable.t24,wvl);
phi1_realtime = 3*phi2_realtime;
cos_phi_2_r = cosd(phi2_realtime); sin_phi_2_r = sind(phi2_realtime);
cos_phi_1_r = cosd(phi1_realtime); sin_phi_1_r = sind(phi1_realtime);
sin_phi1_s2_r = sin_phi_1_r.*sin_phi_2_r;
sin_phi1_c2_r = sin_phi_1_r.*cos_phi_2_r;

%% 生成样品参数和仿真测量光谱
[M_s, N, C, S] = rcwa_film(wvl/1000,0.8,65); % 微米(行向量)，微米，角度
Iout = instruModel(phi1_realtime,phi2_realtime,Iin,M_s,0,num,configuration); % 使用测量时温度

% 由于是随机噪声，重复十次取平均值
for i=1:9
    SNR = 10+5*(i-1);
    Iout_1 = awgn(Iout, SNR, 'measured');

    % 压缩感知方法
    [N_cs, C_cs, S_cs, M_basis, I_cs] = compressed_sensing(Iout_1, cos_phi_1, ...
        sin_phi1_s2, sin_phi1_c2, 5, hyper_paras1(1), hyper_paras1(2), length(wvn));

    % 智能重构算法
    % 存储数据与python交互
    T1 = table(hyper_paras2);
    writetable(T1,"hyper.csv");
    % 角度单位：度°
    T2 = table(Iout_1', phi1', phi2', wvl'*0.001,...
        'VariableNames',{'Y','d1','d2','wvl_um'});
    writetable(T2,"data.csv");
    T3 = table(M_basis);
    writetable(T3,"M_basis.csv");

    % 读取结果
    disp('程序将暂停，使用Net进行求解，请输入“y”以继续运行。');
    user_input = '';
    while ~strcmp(user_input, 'y')
        user_input = input('请输入: ', 's');
    end
    disp('程序继续运行。');

    % 继续执行后续代码
    T_result_Net = readtable("result.csv");
    T_loss_Net = readtable("loss.csv");
    N_Net = (T_result_Net.out_N)';
    C_Net = (T_result_Net.out_C)';
    S_Net= (T_result_Net.out_S)';

    [~,~,~,value_cs(i)] = RMSE(N_cs,N,C_cs,C,S_cs,S,1:length(wvn));
    [~,~,~,value_net(i)] = RMSE(N_Net,N,C_Net,C,S_Net,S,1:length(wvn));
end

%%
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
plot(10:5:50, value_cs,'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(10:5:50, value_net,'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
ylim([0.05,0.27])
leg = legend('CS','Our Method');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 4;

set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig3_10','-dpng','-r600')

