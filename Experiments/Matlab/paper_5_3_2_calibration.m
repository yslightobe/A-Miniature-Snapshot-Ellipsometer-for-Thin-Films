%% 进行仪器的原位同步校准
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
mea_nan = 40;
measeries = mea_nan:num-mea_nan;
meawvl = wvl(measeries);
meanum = length(meawvl);
hyper_paras1 = [0.01, 0.005];
hyper_paras2 = [0.01, 0.005, 0.3, 1, 0];
win = blackman(num)';

ch0 = [344, 358];
ch1 = [362, 372];
ch2 = [373, 387];
ch3 = [388, 402];
ch4 = [403, 417];

%% 读取校准数据
% 0-0-0-0
Iin_1 = readData("..\data\2_Iin_with_4Components", wvl);
% 0-45-0-0
Iout_1 = readData("..\data\4_Iout_straight_air_0-45-0-0", wvl);
% 0-45-0-45
Iout_2 = readData("..\data\5_Iout_straight_air_0-45-0-45", wvl);
% 0-45-0-135
Iout_3 = readData("..\data\6_Iout_straight_air_0-45-0-45+90", wvl);

% ART校准Iin
Iin = (Iout_2 + Iout_3)*2;

A_1 = zero_clean_auto(Iout_1,Iin_1,'blackman',num, ch0);
Iout_1 = Iout_1./Iin_1./A_1;
A_2 = zero_clean_auto(Iout_2,Iin,'blackman',num, ch0);
Iout_2 = Iout_2./Iin./A_2;
A_3 = zero_clean_auto(Iout_3,Iin,'blackman',num, ch0);
Iout_3 = Iout_3./Iin./A_3;

%% 校准相位延迟量
[e_delta1, e_delta1minus2, e_delta1plus2, e_delta2] = channel_calibration_auto(Iout_1, Iout_2, ch0, ch2, ch3, ch4);
% 校准延迟量
d1_cali = unwrap(angle(e_delta1)); % rad
d1_cali = -d1_cali*180/pi;
d2_cali = unwrap(angle(e_delta2)); % rad
d2_cali = -d2_cali*180/pi;
cos_d1 = cosd(d1_cali); sin_d1 = sind(d1_cali);
cos_d2 = cosd(d2_cali); sin_d2 = sind(d2_cali);

% figure(1)
% plot(wvn(measeries), cos_d1(measeries), wvn(measeries), cos_d2(measeries))
% legend("d1", "d2")

%% 提取空气参数
Iout_film = Iout_3;
[N_cs, C_cs, S_cs, M_basis, I_cs] = compressed_sensing_diffConfig("0-45-0-135", ...
    Iout_film(measeries), cos_d1(measeries),sin_d1(measeries).*sin_d2(measeries),...
    sin_d1(measeries).*cos_d2(measeries), 5, hyper_paras1(1), hyper_paras1(2), length(measeries));

% figure(2)
% plot(wvn(measeries), zeros(1,meanum),'--', wvn(measeries), N_cs)
% hold on
% plot(wvn(measeries), ones(1, meanum),'--', wvn(measeries), C_cs)
% plot( wvn(measeries), zeros(1,meanum),'--', wvn(measeries), S_cs)
% legend("N-GT","N-CS","C-GT","C-CS","S-GT","S-CS")

%% 存储数据与python交互
T1 = table(hyper_paras2);
writetable(T1,"hyper.csv");
% 角度单位：度°
T2 = table(Iout_film(measeries)', d1_cali(measeries)', d2_cali(measeries)', wvl(measeries)'*0.001,...
    'VariableNames',{'Y','d1','d2','wvl_um'});
writetable(T2,"data.csv");
T3 = table(M_basis);
writetable(T3,"M_basis.csv");

%% 读取结果
disp('程序将暂停，使用Net进行求解，请输入“y”以继续运行。');
user_input = '';
while ~strcmp(user_input, 'y')
    user_input = input('请输入: ', 's');
end
disp('程序继续运行。');
% 继续执行后续代码
T_result_Net = readtable("result.csv");
T_loss_Net = readtable("loss.csv");
N_Net = T_result_Net.out_N';
C_Net = T_result_Net.out_C';
S_Net = T_result_Net.out_S';
I_Net = T_result_Net.I_predict';
loss_DIP = T_loss_Net.loss_time_eps;

createfigure(wvn(measeries)'*1000,N_Net',zeros(1,meanum)',C_Net',ones(1,meanum)', ...
    S_Net',zeros(1,meanum)',1:meanum,true, [-0.01,0.1], true)

[~,~,~,value_Net] = RMSE(N_Net,zeros(1,meanum),C_Net,ones(1,meanum),S_Net,zeros(1,meanum),1:meanum);

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,1.5*16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(wvn(measeries)*1000, cos_d1(measeries),'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-');
hold on
plot(wvn(measeries)*1000, cos_d2(measeries),'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','--');
xlim([1.25,2.2])
ylim([-1.2,1.6])
leg = legend('$\cos(\varphi_1)$','$\cos(\varphi_2)$','Interpreter','latex');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 2;
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig_5_7_a','-dpng','-r600')


