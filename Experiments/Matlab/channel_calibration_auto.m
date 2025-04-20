function [e_delta1, e_delta1minus2, e_delta1plus2, e_delta2] = channel_calibration_auto(Iout1, Iout2, ch0, ch2, ch3, ch4)
% Iout1: 0-45-0-0的空气光谱
% Iout2: 0-45-0-45的空气光谱
n = length(Iout1); % 计算采样点数
win = blackman(n)'; % 窗函数，备选：hamming hann  blackman barthannwin

f1 = fftshift(fft(Iout1.*win)); % 计算加窗fft(Iin)
f2 = fftshift(fft(Iout2.*win)); % 计算加窗fft(Iin)

figure;
plot(abs(f1), 'LineWidth', 1.6);
title('please choose the 0 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% 提取光强
Ch0_1 = zeros(1, n);
for i = round(ch0(1)):round(ch0(2))
    Ch0_1(i) = f1(i);
end
ch0_1 = ifft(ifftshift(Ch0_1))./win;
% close
% 提取3频延迟量
figure;
plot(abs(f1), 'LineWidth', 1.6);
title('please choose the 3 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
Ch3_1 = zeros(1, n);
for i = round(ch3(1)):round(ch3(2))
    Ch3_1(i) = f1(i);
end
ch3_1 = ifft(ifftshift(Ch3_1))./win;
close;

figure;
plot(abs(f2), 'LineWidth', 1.6);
title('please choose the 0 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
% 提取光强
Ch0_2 = zeros(1, n);
for i = round(ch0(1)):round(ch0(2))
    Ch0_2(i) = f2(i);
end
ch0_2 = ifft(ifftshift(Ch0_2))./win;
% close
% 提取2频延迟量
figure;
plot(abs(f2), 'LineWidth', 1.6);
title('please choose the 2 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
Ch2_2 = zeros(1, n);
for i = round(ch2(1)):round(ch2(2))
    Ch2_2(i) = f2(i);
end
ch2_2 = ifft(ifftshift(Ch2_2))./win;
close;
% 提取4频延迟量
figure;
plot(abs(f2), 'LineWidth', 1.6);
title('please choose the 4 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
Ch4_2 = zeros(1, n);
for i = round(ch4(1)):round(ch4(2))
    Ch4_2(i) = f2(i);
end
ch4_2 = ifft(ifftshift(Ch4_2))./win;
close;

e_delta1 = 2*ch3_1./ch0_1;
e_delta1minus2 = 4*ch2_2./ch0_2;
e_delta1plus2 = -4*ch4_2./ch0_2;
% e_delta2 = e_delta1plus2./e_delta1;
e_delta2 = (e_delta1plus2./e_delta1 + e_delta1./e_delta1minus2)/2;
end