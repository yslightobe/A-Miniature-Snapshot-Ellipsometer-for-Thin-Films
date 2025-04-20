% 读取数据文件
filename = '../data/薄膜厚度.txt'; % 文件路径

% 使用readtable读取数据（支持更多格式控制）
try
    % 读取空格分隔的文本文件，无标题行，自动填充NaN
    data_table = readtable(filename, ...
        'Delimiter', '\t', ...       % 指定分隔符为空格
        'ReadVariableNames', false,... % 无列标题
        'NumHeaderLines', 0);       % 跳过标题行
    
    % 将表格转换为数值矩阵
    raw_data = table2array(data_table);
catch ME
    error('文件读取失败！请检查文件路径和格式。错误信息：%s', ME.message);
end

% 初始化结果数组
num_samples = size(raw_data, 1); % 样品数量（行数）
avg_values = zeros(num_samples, 1);
three_sigma = zeros(num_samples, 1);
uncertainty = zeros(num_samples, 1);
valid_counts = zeros(num_samples, 1); % 记录有效测量次数

% 循环处理每个样品（每行）
for i = 1:num_samples
    measurements = raw_data(i, :); % 当前行的原始数据
    
    % 移除NaN值（处理因列数不一致导致的NaN）
    clean_data = measurements(~isnan(measurements));
    
    n = length(clean_data); % 有效测量次数
    
    % 如果有效数据不足（如n=0），标记为NaN
    if n < 1
        avg_values(i) = NaN;
        three_sigma(i) = NaN;
        uncertainty(i) = NaN;
    else
        avg_values(i) = mean(clean_data); % 平均值
        std_dev = std(clean_data); % 标准差
        three_sigma(i) = 3 * std_dev; % 3σ
        uncertainty(i) = std_dev / sqrt(n); % 不确定度（标准误差）
    end
    
    valid_counts(i) = n; % 记录有效测量次数
end

% 创建结果表格
results = table(...
    (1:num_samples)', ... % 样品编号（1到num_samples）
    avg_values, ...
    three_sigma, ...
    uncertainty, ...
    valid_counts, ...
    'VariableNames', {'Sample', 'Average', '3σ', 'Uncertainty', 'ValidMeasurements'});

% 显示结果
disp('薄膜厚度测量结果统计（处理NaN后）：');
disp(results);

% 保存结果到文件（可选）

% writetable(results, '../data/thickness_results_processed.txt');
% disp('结果已保存为 thickness_results_processed.txt');

