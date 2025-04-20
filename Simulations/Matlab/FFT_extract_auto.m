function [N_extracted,C_extracted,S_extracted] = FFT_extract_auto(Iin,Iout,wvn,ch0,ch1,ch2,ch3,e_delta1,e_delta1minusdelta2)
%  归一化测量光谱

inten = Iout./Iin; % 消除背景光谱
% max_inten = max(inten); % 光谱数据最大值
% inten = inten/max_inten; % 归一化测量光谱
n = length(wvn); % 计算采样点数
win = blackman(n)'; % 窗函数，备选：hamming hann  blackman barthannwin
Inten = fftshift(fft(inten.*win)); % 计算加窗fft(Iin)
sample_points = length(wvn);

% 画出校准光谱的OPD域图并选频
% figure;
% plot(abs(Inten), 'LineWidth', 1.6);
% title('please choose the 0 channel')
% set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% 提取光强
% [f_cat_f, ~] = ginput(2);
Ch0_f = zeros(1, sample_points);
for i = ch0(1):ch0(2)
    Ch0_f(i) = Inten(i);
end
ch0_f = ifft(ifftshift(Ch0_f))./win;
% close
% 提取1频延迟量
% figure;
% plot(abs(Inten), 'LineWidth', 1.6);
% title('please choose the 1 channel')
% set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% [f_cat_f, ~] = ginput(2);
Ch1_f = zeros(1, sample_points);
for i = ch1(1):ch1(2)
    Ch1_f(i) = Inten(i);
end
ch1_f = ifft(ifftshift(Ch1_f))./win;
% close;

% 提取2频延迟量
% figure;
% plot(abs(Inten), 'LineWidth', 1.6);
% title('please choose the 2 channel')
% set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% [f_cat_f, ~] = ginput(2);
Ch2_f = zeros(1, sample_points);
for i = ch2(1):ch2(2)
    Ch2_f(i) = Inten(i);
end
ch2_f = ifft(ifftshift(Ch2_f))./win;
% close;

% 提取3频延迟量
% figure;
% plot(abs(Inten), 'LineWidth', 1.6);
% title('please choose the 4 channel')
% set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% [f_cat_f, ~] = ginput(2);
Ch3_f = zeros(1, sample_points);
for i = ch3(1):ch3(2)
    Ch3_f(i) = Inten(i);
end
ch3_f = ifft(ifftshift(Ch3_f))./win;
% close;

N_extracted = (-real(2*ch3_f./ch0_f./e_delta1))';
C_extracted = (real(-4*ch2_f./ch0_f./e_delta1minusdelta2))';
S_extracted = (imag(-4*ch2_f./ch0_f./e_delta1minusdelta2))';
end