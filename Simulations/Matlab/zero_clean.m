% 用于光谱零频的清洗
% created by：杨世龙
% date：2023.6
function A = zero_clean(Iout, Iin, window_name, sample_points)
% Iout, Iin: 光谱数据
eval(strcat('win = ', window_name, "(sample_points);"));
f = fftshift(fft(8*Iout./Iin.*win));
figure;
plot(abs(f), 'LineWidth', 1.6);
title('please choose the 0 channel')
set(gca, 'Fontsize', 16, 'Fontname', 'Times New Roman')
[f_cat, ~] = ginput(2);
ch0 = zeros(sample_points,1);
for i = round(f_cat(1)):round(f_cat(2))
    ch0(i) = f(i);
end
temp = real(ifft(ifftshift(ch0))./win);
close
A = temp/2;
end