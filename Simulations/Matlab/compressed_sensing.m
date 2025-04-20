% 用于NCS的压缩感知(2017文献)解调
% created by：YSL
% date：2023.6
function  [N, C, S, M_basis, Igenerate] = compressed_sensing(Y, cos_delta_1,sin_delta1_s2,sin_delta1_c2, L, DCT, beta, sample_points)
% 用于通道型Stokes偏振仪的计算
% Iout, Iin：测量光谱，列向量（N×1）
% N, C, S：薄膜参数，列向量（N×1）
% delta_2, delta_1Minus2, delta_1Plus2：仪器里的波片延迟量
% L：勒让德多项式阶数，正整数
% DCT：DCT矩阵频率限，0-1小数
% beta：损失函数中稀疏度的比例
% sample_points：采样点数，正整数
%% 压缩感知方法求解stokes参量
% 构造三角函数矩阵
M_tri_N = diag(cos_delta_1);
M_tri_C = diag(sin_delta1_s2);
M_tri_S = diag(sin_delta1_c2);

% 构造DCT矩阵
M_dct = zeros(sample_points, sample_points);
f_set = round((L+sample_points)*DCT); % 设置最高频率，防止高频噪声的扰动
for i = 1:sample_points
    for j = 1:sample_points
        if i==1
            M_dct(1, :) = (1/sample_points)^0.5;
        else
            M_dct(i, j) = (2/sample_points)^0.5*cos(pi/(2*sample_points)*(i-1)*(2*j-1));
        end
    end
end
% 构造Legendre矩阵
M_legendre = zeros(sample_points, L);
for i = 1:L
    M_temp = legendre(i-1, linspace(-1,1,sample_points));
    M_legendre(:, i) = M_temp(1, :)';
end
% 构造support矩阵
M_basis = [M_legendre M_dct];

M_final = [M_tri_N*M_basis M_tri_C*M_basis M_tri_S*M_basis];
zero_set = zeros([L+sample_points-f_set,1]);
% 对此压缩感知问题进行求解
cvx_clear
cvx_begin
    variable x_optim_0(f_set)
    variable x_optim_1(f_set)
    variable x_optim_2(f_set)
    Igenerate = 0.25*(M_final*[x_optim_0;zero_set;x_optim_1;zero_set;x_optim_2;zero_set]+1);
    minimize(beta*(norm(x_optim_0, 1)+norm(x_optim_1, 1)+norm(x_optim_2, 1))...
        + norm(Igenerate-Y', 2));
cvx_end
y_out = blkdiag(M_basis, M_basis, M_basis)*[x_optim_0;zero_set;x_optim_1;zero_set;x_optim_2;zero_set];
%% 输出结果
Igenerate = Igenerate';
N = y_out(1:sample_points, 1)';
C = y_out(sample_points+1:2*sample_points, 1)';
S = y_out(2*sample_points+1:3*sample_points, 1)';
end