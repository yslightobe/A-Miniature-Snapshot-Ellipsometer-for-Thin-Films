% 用于光谱零频的清洗
% created by：杨世龙
% date：2023.6
function A = zero_clean_auto(Iout, Iin, window_name, sample_points, ch0)
% Iout, Iin: 光谱数据
eval(strcat('win = ', window_name, "(sample_points);"));
f = fftshift(fft(Iout./Iin.*win'));
figure;
plot(abs(f), 'LineWidth', 1.6);
title('please choose the 0 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
Ch0 = zeros(1,sample_points);
for i = round(ch0(1)):round(ch0(2))
    Ch0(i) = f(i);
end
temp = real(ifft(ifftshift(Ch0))./win');
close
A = temp;
end