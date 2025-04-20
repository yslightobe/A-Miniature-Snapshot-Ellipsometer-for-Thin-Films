function [spectra, spec_mat_raw, wvl_raw] = read_data(file_root, wvl)
% created by：YSL
% date：2024.3
% 读取光谱信息，并整形
% file_root: data path
% wvl: wavelength vector
% spectra: average data after fit
% spec_mat_raw: matrix consist of raw spectrum series
% wvl_raw: raw wavelength from spectrometer
list = dir(strcat(file_root, "\*.csv"));
num = length(list);
spec_mat_raw = zeros(4011, num);
for i=1:num
    T1 = readtable(strcat(file_root, '\', list(i).name));
    spec_mat_raw(:,i) = T1.Var2;
end
wvl_raw = T1.Var1;
spectra_raw = mean(spec_mat_raw, 2);
spectra = spline(wvl_raw, spectra_raw, wvl);
x = ['***共采集', num2str(num), '组光谱数据***'];
disp(x)
end