function [spectra, spec_mat_raw, wvl_raw] = readData_2(file_root, wvl)
% code by YSL (modified for natural order with specific filename format)
% time: 202403
% 读取光谱信息，并整形
% file_root: data path
% wvl: wavelength vector
% spectra: average data after fit
% spec_mat_raw: matrix consist of raw spectrum series
% wvl_raw: raw wavelength from spectrometer

% 获取所有CSV文件列表
list = dir(fullfile(file_root, "Spectrum_*.csv"));  % 根据文件名前缀筛选
num = length(list);

% 提取时间戳部分并排序（自然顺序）
numbers = cell(num, 1);
for i = 1:num
    filename = list(i).name;

    pattern = '_([0-9]+\.[0-9]+) ms\.csv$';
    matches = regexpi(filename, pattern, 'tokens');
    
    if ~isempty(matches{1})
        time_str = matches{1}{1};
        numbers{i} = str2double(time_str);
    else
        error('文件名中未找到有效时间戳，请检查文件命名格式！');
    end
end

% 将数字转换为数值数组并排序
[~, idx] = sort(cell2mat(numbers));
list = list(idx);  % 重新排列文件列表顺序

% 初始化矩阵（需先读取第一个文件获取维度）
first_file = fullfile(file_root, list(1).name);
T1 = readtable(first_file, ...
    "NumHeaderLines", 26, ...
    "VariableNamingRule", "preserve", ...
    "VariableNamesLine", 26);
spec_mat_raw = zeros(size(T1.Intensity, 1), num);

% 读取所有光谱数据
for i = 1:num
    current_file = fullfile(file_root, list(i).name);
    T1 = readtable(current_file, ...
        "NumHeaderLines", 26, ...
        "VariableNamingRule", "preserve", ...
        "VariableNamesLine", 26);
    spec_mat_raw(:, i) = T1.Intensity;
end

% 获取波长数据（假设所有文件的波长相同）
wvl_raw = T1.("Wavelength(nm)");

% 计算平均光谱并插值到指定波长
spectra_raw = mean(spec_mat_raw, 2);
spectra = spline(wvl_raw', spectra_raw', wvl);

% 显示信息
disp(strcat('***共采集 ', num2str(num), ' 组光谱数据（按时间戳排序）***'));
end
