%% 3.2.2 基矩阵示意图
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
num = 30;

%% 构造基矩阵
% 构造DCT矩阵
M_dct = zeros(num, num);
for i = 1:num
    for j = 1:num
        if i==1
            M_dct(1, :) = (1/num)^0.5;
        else
            M_dct(i, j) = (2/num)^0.5*cos(pi/(2*num)*(i-1)*(2*j-1));
        end
    end
end
% 构造Legendre矩阵
M_legendre = zeros(num, 5);
for i = 1:5
    M_temp = legendre(i-1, linspace(-1,1,num));
    M_legendre(:, i) = M_temp(1, :)';
end
% 构造support矩阵
M_basis = [M_legendre M_dct];

%% 绘制图像
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
imagesc(M_basis)
cb = colorbar;
cb.TickLabels = [];
set(cb, 'LineWidth', 1)
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
% print('fig3_2_a','-dpng','-r600')
