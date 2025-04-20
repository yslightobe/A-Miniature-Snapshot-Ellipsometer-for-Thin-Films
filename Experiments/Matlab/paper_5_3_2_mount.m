%% 装调过程中各频道强度的变化
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
hyper_paras1 = [0.02, 0.02];
hyper_paras2 = [0.02, 0.02, 0.3];
win = blackman(num)';

ch0 = [344, 358];
ch1 = [362, 372];
ch2 = [373, 387];
ch3 = [388, 402];
ch4 = [403, 417];

%% 读取数据
% 光强原始响应
Iin_0 = readData("D:\TianyiPanFiles\PostGraduateFiles\Documents" + ...
    "\实验室文献\硕士毕业论文\code\experiment\data\1_SourceSpectrum", wvl);

% fig5-5(a)
[~, spec_mat, wvl_raw] = readData_2("D:\TianyiPanFiles\PostGraduateFiles\" + ...
    "Documents\实验室文献\硕士毕业论文\code\experiment\data\realtime1_0-0-null-0", wvl);

ch3_r1_null = zeros(1, size(spec_mat,2));
for i=1:size(spec_mat,2)
    spectrum = (spline(wvl_raw', spec_mat(:,i)', wvl));
    f = abs(fftshift(fft(spectrum.*win)));
%     figure(1);
%     plot(OPD, f);
    ch3_r1_null(i) = max(f(ch3(1):ch3(2)));
end
% figure(2)
% plot(ch3_r1)

% fig5-5(b)
[~, spec_mat, wvl_raw] = readData_2("D:\TianyiPanFiles\PostGraduateFiles\" + ...
    "Documents\实验室文献\硕士毕业论文\code\experiment\data\realtime2_0-0-0-0", wvl);

ch1_r2 = zeros(1, size(spec_mat,2));
for i=1:size(spec_mat,2)
    spectrum = (spline(wvl_raw', spec_mat(:,i)', wvl));
    f = abs(fftshift(fft(spectrum.*win)));
%     figure(3);
%     plot(OPD, f);
    ch1_r2(i) = max(f(ch1(1):ch1(2)));
end
% figure(4)
% plot(ch1_r2)

%%
figure(1)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(ch3_r1_null/max(ch3_r1_null),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-3,190])
ylim([-0.05,1.1])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
print('fig_5_5_a','-dpng','-r600')

figure(2)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(ch1_r2/max(ch1_r2),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-3,400])
ylim([-0.05,1.1])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
print('fig_5_5_b','-dpng','-r600')

%%

% fig5-6(a)
[~, spec_mat, wvl_raw] = readData_2("D:\TianyiPanFiles\PostGraduateFiles\" + ...
    "Documents\实验室文献\硕士毕业论文\code\experiment\data\realtime3_0-45-0-0", wvl);

ch3_r1 = zeros(1, size(spec_mat,2));
for i=1:size(spec_mat,2)
    spectrum = (spline(wvl_raw', spec_mat(:,i)', wvl));
    f = abs(fftshift(fft(spectrum.*win)));
%     figure(5);
%     plot(OPD, f);
    ch3_r1(i) = max(f(ch3(1):ch3(2)));
end
% figure(6)
% plot(ch3_r1)

% fig5-6(b)
[~, spec_mat, wvl_raw] = readData_2("D:\TianyiPanFiles\PostGraduateFiles\" + ...
    "Documents\实验室文献\硕士毕业论文\code\experiment\data\realtime4_0-45-0-45", wvl);

ch2_a = zeros(1, size(spec_mat,2));
ch3_a = zeros(1, size(spec_mat,2));
ch4_a = zeros(1, size(spec_mat,2));
for i=1:size(spec_mat,2)
    spectrum = (spline(wvl_raw', spec_mat(:,i)', wvl));
    f = abs(fftshift(fft(spectrum.*win)));
%     figure(7);
%     plot(OPD, f);
    ch2_a(i) = max(f(ch2(1):ch2(2)));
    ch3_a(i) = max(f(ch3(1):ch3(2)));
    ch4_a(i) = max(f(ch4(1):ch4(2)));
end
% figure(8)
% plot(1:size(spec_mat,2),ch2_a,1:size(spec_mat,2),ch3_a,1:size(spec_mat,2),ch4_a)
% legend("2频","3频","4频")

%%
figure(3)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(ch3_r1/max(ch3_r1),'LineWidth',1.2,'Color',[0 0 0],'LineStyle','-');
xlim([-3,220])
ylim([-0.05,1.1])
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
print('fig_5_6_a','-dpng','-r600')

figure(4)
set(gcf,'unit', 'centimeters', 'position', [25,9.5,16/2.5,9/2.5]);
% 使用 plot 的矩阵输入创建多行
plot(ch2_a/max(ch2_a),'LineWidth',1.2,'Color',[0 0.447 0.741],'LineStyle','-'); hold on
plot(ch3_a/max(ch3_a),'LineWidth',1.2,'Color',[0.851 0.325 0.098],'LineStyle','-');
plot(ch4_a/max(ch4_a),'LineWidth',1.2,'Color',[0.467 0.675 0.188],'LineStyle','-');
xlim([-3,623])
ylim([-0.05,1.1])
leg = legend('±2','±3','±4');
leg.ItemTokenSize = [12,5]; leg.Box = 'off'; leg.NumColumns = 3; leg.Location = 'east';
set(gca,'LooseInset',[0.01 0.01 0.01 0.01],"LineWidth",1,"Fontname","Times New Roman", ...
    "Fontsize",12,"XMinorTick","on","YMinorTick","on",xticklabels=[],yticklabels=[])
print('fig_5_6_b','-dpng','-r600')