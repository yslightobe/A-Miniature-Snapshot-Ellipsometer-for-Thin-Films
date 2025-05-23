import torch
import torch.optim
import torch.nn as nn
import numpy as np
import pandas as pd
import random
import time
from memory_profiler import profile
import os
from torch import dtype


def read_col(file_path, col):
    """用于读取csv文件的数据
    :param file_path: 文件地址（字符串）
    :param col: 欲读取数据的列名（字符串）
    :return: 返回行向量数据
    """
    data = pd.read_csv(file_path)
    coldata = data[col]
    coldata = np.array(coldata, dtype=np.float64)
    outdata = coldata.reshape((1, coldata.size))
    return outdata


def setup_seed(seed):
    torch.manual_seed(seed)
    # torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    # torch.backends.cudnn.deterministic = True


def RetardanceT(wvl, T, d=1.5):
    """波片延迟量温度特性模型

    :param wvl: 波长，单位：μm
    :param d: 波片厚度，单位：mm
    :param T: 环境温度，单位：℃
    """
    ne = torch.sqrt(1 + 0.665721 * wvl ** 2 / (wvl ** 2 - 0.06 ** 2) + 0.503511 * wvl ** 2 / (
            wvl ** 2 - 0.106 ** 2) + 0.214792 * wvl ** 2 / (wvl ** 2 - 0.119 ** 2) + 0.539173 * wvl ** 2 / (
                            wvl ** 2 - 8.792 ** 2) + 1.8076613 * wvl ** 2 / (wvl ** 2 - 19.7 ** 2))
    no = torch.sqrt(1 + 0.663044 * wvl ** 2 / (wvl ** 2 - 0.06 ** 2) + 0.517852 * wvl ** 2 / (
            wvl ** 2 - 0.106 ** 2) + 0.175912 * wvl ** 2 / (wvl ** 2 - 0.119 ** 2) + 0.56538 * wvl ** 2 / (
                            wvl ** 2 - 8.844 ** 2) + 1.675299 * wvl ** 2 / (wvl ** 2 - 20.742 ** 2))
    lambda_ig_o = 1240 / 10.30 * 10 ** (-3)
    R_o = wvl ** 2 / (wvl ** 2 - lambda_ig_o ** 2)
    dno = (-61.1840 * R_o + 43.999 * R_o ** 2) / (2 * no) * 10 ** (-6)
    lambda_ig_e = 1240 / 10.30 * 10 ** (-3)
    R_e = wvl ** 2 / (wvl ** 2 - lambda_ig_e ** 2)
    dne = (-70.1182 * R_e + 49.2875 * R_e ** 2) / (2 * ne) * 10 ** (-6)
    delta0 = 2 * np.pi * d / wvl * (ne - no) * 10 ** 3
    K = 13.37 * 10 ** (-6) + (dne - dno) / (ne - no)
    d_delta = delta0 * (T - 25) * K
    return d_delta


def polarizer(theta):
    """用于生成偏振片穆勒矩阵的函数
    :param theta: 偏振片的方位角（弧度）
    """
    M_p = torch.tensor([[0.5, 0.5*torch.cos(2*theta), 0.5*torch.sin(2*theta), 0],
                        [0.5*torch.cos(2*theta), 0.5*torch.cos(2*theta)**2, 0.5*torch.sin(2*theta)*torch.cos(2*theta), 0],
                        [0.5*torch.sin(2*theta), 0.5*torch.sin(2*theta)*torch.cos(2*theta), 0.5*torch.sin(2*theta)**2, 0],
                        [0, 0, 0, 0]], dtype=torch.float64)
    return M_p


def waveplate(theta, delta):
    """用于生成波片穆勒矩阵的函数
    :param theta: 波片的方位角（弧度），尺寸为[1, 1, 1]
    :param delta: 波片的延迟量（弧度），尺寸为[1, 1, n]
    """
    # 确保theta和delta的形状正确
    theta = theta.squeeze(0).squeeze(0)  # 将theta从[1, 1, 1]变为标量
    delta = delta.squeeze(0).squeeze(0)  # 将delta从[1, 1, n]变为[n]

    # 预计算常用的三角函数
    cos_delta = torch.cos(delta)
    sin_delta = torch.sin(delta)
    sin_2theta = torch.sin(2 * theta)
    cos_2theta = torch.cos(2 * theta)
    cos_2delta = torch.cos(2 * delta)

    # 生成穆勒矩阵的各个元素
    M_R_00 = torch.ones_like(delta)
    M_R_11 = 1 - (1 - cos_delta) * sin_2theta**2
    M_R_12 = (1 - cos_2delta) * sin_2theta * cos_2theta
    M_R_13 = -sin_delta * sin_2theta
    M_R_21 = (1 - cos_delta) * sin_2theta * cos_2theta
    M_R_22 = 1 - (1 - cos_delta) * cos_2theta**2
    M_R_23 = sin_delta * cos_2theta
    M_R_31 = sin_delta * sin_2theta
    M_R_32 = -sin_delta * cos_2theta
    M_R_33 = cos_delta

    # 样品矩阵初始化
    M_R = torch.zeros(len(delta), 4, 4, dtype=torch.float64)

    # 填充穆勒矩阵的各个元素
    M_R[:, 0, 0] = M_R_00
    M_R[:,1, 1] = M_R_11
    M_R[:, 1, 2] = M_R_12
    M_R[:, 1, 3] = M_R_13
    M_R[:, 2, 1] = M_R_21
    M_R[:, 2, 2] = M_R_22
    M_R[:, 2, 3] = M_R_23
    M_R[:, 3, 1] = M_R_31
    M_R[:, 3, 2] = M_R_32
    M_R[:, 3, 3] = M_R_33

    return M_R


def Sample(N, C, S):
    """用于生成样品矩阵的函数
    :param N: 样品的偏振参数1
    :param C: 样品的偏振参数2
    :param S: 样品的偏振参数3
    """

    # 样品矩阵初始化
    M_s = torch.zeros(N.shape[2], 4, 4, dtype=torch.float64)
    M_s_ones = torch.ones_like(N[0, 0, :])

    # 填充样品矩阵的各个元素
    M_s[:, 0, 0] = M_s_ones
    M_s[:, 0, 1] = N[0, 0, :]
    M_s[:, 1, 0] = N[0, 0, :]
    M_s[:, 1, 1] = M_s_ones
    M_s[:, 2, 2] = C[0, 0, :]
    M_s[:, 2, 3] = S[0, 0, :]
    M_s[:, 3, 2] = -S[0, 0, :]
    M_s[:, 3, 3] = C[0, 0, :]

    return M_s


class UNN(nn.Module):
    def __init__(self, in_ch, hidden_ch, out_ch):
        super(UNN, self).__init__()
        # 输出基系数
        self.StackLayer = nn.Sequential(
            nn.Linear(in_ch, hidden_ch, dtype=torch.float64),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_ch, hidden_ch, dtype=torch.float64),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_ch, out_ch, dtype=torch.float64)
        )
        # 输出方位角误差
        self.layer1 = nn.Sequential(
            nn.Tanh(),
            nn.Linear(out_ch, 4, dtype=torch.float64)
        )
        # 输出波片温度控制系数
        self.layer2 = nn.Sequential(
            nn.Tanh(),
            nn.Linear(4, 1, dtype=torch.float64)
        )

    def forward(self, x):
        out1 = self.StackLayer(x)
        out2 = self.layer1(out1) * np.pi / 180
        out3 = self.layer2(out2)
        return out1, out2, out3


@profile  # 添加此装饰器
def SFCNet(wvl, M_basis, Y, d1, d2, num=2, num2=400,
           f_set=0.05, beta=0.01, gama=1, seed=5, range1=0.05, range2=0.95, align=False, temperature=False):
    # 不分组
    series_num = int(M_basis.shape[2] * f_set)
    zero_col = torch.zeros([1, 1, M_basis.shape[2] - series_num])
    sample_points = int(Y.shape[2])
    setup_seed(seed)
    myNet = UNN(sample_points, 50, 3 * series_num)
    start_lr = 0.001
    params = list(myNet.StackLayer.parameters()) + list(myNet.layer1.parameters())
    optimizer = torch.optim.Adam(params, lr=start_lr)
    running_loss = np.zeros(num * num2)
    Iin = torch.tensor([[[1],[0],[0],[0]]], dtype=torch.float64)
    # 调整不同的理论角度
    angle_R1 = torch.tensor([[[np.pi/4]]]) # np.pi/4 3*np.pi/4
    angle_R2 = torch.tensor([[[0]]]) # 0 np.pi/2
    angle_A = torch.tensor([[[np.pi/4]]]) # np.pi/4 3*np.pi/4
    start = time.time()

    for epoch in range(num):
        for iter_num in range(num2):
            optimizer = optimizer
            optimizer.zero_grad()
            output1, output2, output3 = myNet(Y)

            # 预测可能存在的方位角误差（P-R1-R2-A）
            if align:
                beta1 = output2[:1, :1, 0]
                alpha1 = output2[:1, :1, 1]
                beta2 = output2[:1, :1, 2]
                alpha2 = output2[:1, :1, 3]
            else:
                beta1 = torch.tensor(0)
                alpha1 = torch.tensor(0)
                beta2 = torch.tensor(0)
                alpha2 = torch.tensor(0)


            # 预测波片温度控制系数
            if temperature:
                T = (output3[:1, :1, 0] + 1) * 20 - 20 + 24
                d_delta2 = RetardanceT(wvl, T) * 10
                d_delta1 = 3 * d_delta2
            else:
                d_delta1 = 0
                d_delta2 = 0
                T = torch.tensor(24)

            temp1 = torch.cat([output1[:1, :1, 0:series_num], zero_col], dim=2)
            temp2 = torch.cat([output1[:1, :1, series_num:2 * series_num], zero_col], dim=2)
            temp3 = torch.cat([output1[:1, :1, 2 * series_num:3 * series_num], zero_col], dim=2)

            N_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp1])
            C_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp2])
            S_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp3])

            M_s = Sample(N_g, C_g, S_g)
            M_P = polarizer(beta1)
            M_R1 = waveplate(alpha1+angle_R1, d1+d_delta1)
            M_R2 = waveplate(beta2+angle_R2, d2+d_delta2)
            M_A = polarizer(alpha2+angle_A)

            I_g_temp = M_A @ M_R2 @ M_s @ M_R1 @ M_P @ Iin
            I_g_extract = I_g_temp[:, 0, 0]
            I_g = I_g_extract.unsqueeze(0).unsqueeze(0)

            error = (Y[:1, :1, 1:-1] - I_g[:1, :1, 1:-1])

            loss = (torch.norm(error[:1, :1, int(range1 * sample_points):int(range2 * sample_points)], p=2,
                               dim=2).squeeze(0)
                    + beta * torch.norm(output1, p=1, dim=2).squeeze(0)
                    + gama * torch.norm(N_g ** 2 + C_g ** 2 +
                                        S_g ** 2 - 1, p=2, dim=2).squeeze(0))

            loss.backward()
            optimizer.step()
            running_loss[epoch * num2 + iter_num] = loss.item()

    end = time.time()
    time_cost = end - start
    print('iter num: %d, loss=%f' % (epoch * num2 + iter_num + 1, loss.item()), '\n')
    print('time cost: ', time_cost, 's \n')

    return N_g, C_g, S_g, I_g, beta1, alpha1, beta2, alpha2, running_loss, time_cost, T


@profile  # 添加此装饰器
def SFCNet_GPU(wvl, M_basis, Y, d1, d2, num=2, num2=400,
           f_set=0.05, beta=0.01, gama=1, seed=5, range1=0.05, range2=0.95, align=False, temperature=False, device=None):
    # 这里的实现非常粗糙，只是把张量挪移到GPU上，可以优化的地方还有很多
    # 确保所有张量在GPU上
    wvl = wvl.to(device)
    M_basis = M_basis.to(device)
    Y = Y.to(device)
    d1 = d1.to(device)
    d2 = d2.to(device)

    series_num = int(M_basis.shape[2] * f_set)
    zero_col = torch.zeros([1, 1, M_basis.shape[2] - series_num], device=device)
    sample_points = int(Y.shape[2])
    setup_seed(seed)
    myNet = UNN(sample_points, 50, 3 * series_num).to(device)  # 模型移到GPU
    start_lr = 0.001
    params = list(myNet.StackLayer.parameters()) + list(myNet.layer1.parameters())
    optimizer = torch.optim.Adam(params, lr=start_lr)
    running_loss = torch.zeros(num * num2, device=device)  # 使用GPU张量
    Iin = torch.tensor([[[1],[0],[0],[0]]], dtype=torch.float64, device=device)
    # 调整角度张量到GPU
    angle_R1 = torch.tensor([[[np.pi/4]]], device=device)
    angle_R2 = torch.tensor([[[0.0]]], device=device)
    angle_A = torch.tensor([[[np.pi/4]]], device=device)
    start = time.time()

    for epoch in range(num):
        for iter_num in range(num2):
            optimizer.zero_grad()
            output1, output2, output3 = myNet(Y)

            if align:
                beta1 = output2[:1, :1, 0]
                alpha1 = output2[:1, :1, 1]
                beta2 = output2[:1, :1, 2]
                alpha2 = output2[:1, :1, 3]
            else:
                beta1 = torch.tensor(0.0, device=device)
                alpha1 = torch.tensor(0.0, device=device)
                beta2 = torch.tensor(0.0, device=device)
                alpha2 = torch.tensor(0.0, device=device)

            if temperature:
                T = (output3[:1, :1, 0] + 1) * 20 - 20 + 24
                d_delta2 = RetardanceT(wvl, T) * 10
                d_delta1 = 3 * d_delta2
            else:
                d_delta1 = torch.tensor(0.0, device=device)
                d_delta2 = torch.tensor(0.0, device=device)
                T = torch.tensor(24.0, device=device)

            # 使用GPU张量操作
            temp1 = torch.cat([output1[:1, :1, 0:series_num], zero_col], dim=2)
            temp2 = torch.cat([output1[:1, :1, series_num:2*series_num], zero_col], dim=2)
            temp3 = torch.cat([output1[:1, :1, 2*series_num:3*series_num], zero_col], dim=2)

            N_g = torch.einsum('ijk,ilk->ilj', M_basis, temp1)
            C_g = torch.einsum('ijk,ilk->ilj', M_basis, temp2)
            S_g = torch.einsum('ijk,ilk->ilj', M_basis, temp3)

            M_s = Sample(N_g, C_g, S_g).to(device)
            M_P = polarizer(beta1).to(device)
            M_R1 = waveplate(alpha1+angle_R1, d1+d_delta1).to(device)
            M_R2 = waveplate(beta2+angle_R2, d2+d_delta2).to(device)
            M_A = polarizer(alpha2+angle_A).to(device)

            I_g_temp = M_A @ M_R2 @ M_s @ M_R1 @ M_P @ Iin
            I_g = I_g_temp[:, 0, 0].unsqueeze(0).unsqueeze(0)

            error = Y[:1, :1, 1:-1] - I_g[:1, :1, 1:-1]

            # 使用GPU张量计算
            loss = (
                torch.norm(error[:1, :1, int(range1*sample_points):int(range2*sample_points)], p=2, dim=2).squeeze(0) +
                beta * torch.norm(output1, p=1, dim=2).squeeze(0) +
                gama * torch.norm(N_g**2 + C_g**2 + S_g**2 - 1, p=2, dim=2).squeeze(0)
            )

            loss.backward()
            optimizer.step()
            running_loss[epoch*num2 + iter_num] = loss.item()

    end = time.time()
    time_cost = end - start
    print('iter num: %d, loss=%f' % (epoch * num2 + iter_num + 1, loss.item()), '\n')
    print('time cost: ', time_cost, 's \n')

    return N_g, C_g, S_g, I_g, beta1, alpha1, beta2, alpha2, running_loss.cpu().numpy(), time_cost, T


class Residual(nn.Module):
    def __init__(self, input_channels, num_channels, dowmsampling=False, strides=1):
        super().__init__()
        self.conv1 = nn.Conv1d(input_channels, num_channels,
                               kernel_size=3, padding=1, stride=strides, dtype=torch.float64)
        self.conv2 = nn.Conv1d(num_channels, num_channels,
                               kernel_size=3, padding=1, dtype=torch.float64)
        if dowmsampling:
            self.conv3 = nn.Conv1d(input_channels, num_channels,
                                   kernel_size=3, padding=1, stride=strides, dtype=torch.float64)
        else:
            self.conv3 = None
        self.relu = nn.ReLU()

    def forward(self, X):
        Y = self.relu(self.conv1(X))
        Y = self.conv2(Y)
        if self.conv3:
            X = self.conv3(X)
        Y += X
        return self.relu(Y)

    pass


class ResNet(nn.Module):
    def __init__(self, out_ch):
        super(ResNet, self).__init__()
        b1 = nn.Sequential(nn.Conv1d(1, 64, kernel_size=5, stride=1, padding=2, dtype=torch.float64),
                           nn.ReLU(),
                           nn.MaxPool1d(kernel_size=2, stride=2))

        b2 = nn.Sequential(*self.__resnet_block(64, 64, 3, Fisrtblock=True))
        b3 = nn.Sequential(*self.__resnet_block(64, 128, 4))
        b4 = nn.Sequential(*self.__resnet_block(128, 256, 6))
        b5 = nn.Sequential(*self.__resnet_block(256, 256, 3))
        # self.conv = nn.Sequential(b1, b2, b3, b4, b5, nn.AvgPool1d(2), nn.Flatten())
        self.conv = nn.Sequential(b1, b2, nn.Flatten(), nn.ReLU())
        # self.fc = nn.Linear(3072, out_ch, dtype=torch.float64)
        self.fc = nn.Linear(22400, out_ch, dtype=torch.float64)

        pass

    def forward(self, x):
        x = self.conv(x)
        x = self.fc(x)
        return x.unsqueeze(0)
        pass

    def __resnet_block(self, input_channels, num_channels, num_residuals, Fisrtblock=False):
        blk = []
        for i in range(num_residuals):
            if i == 0 and not Fisrtblock:
                blk.append(Residual(input_channels, num_channels, dowmsampling=True, strides=2))
            else:
                blk.append(Residual(num_channels, num_channels))
        return blk


def SFCNet_Res(wvl, M_basis, Y, d1, d2, num=2, num2=400,
               f_set=0.05, beta=0.01, gama=1, seed=5, range1=0.05, range2=0.95, align=False, temperature=False):
    # 不分组
    series_num = int(M_basis.shape[2] * f_set)
    zero_col = torch.zeros([1, 1, M_basis.shape[2] - series_num])
    sample_points = int(Y.shape[2])
    setup_seed(seed)
    myNet = ResNet(out_ch=3 * series_num + 5)
    start_lr = 0.001
    params = myNet.parameters()
    optimizer = torch.optim.Adam(params, lr=start_lr)
    running_loss = np.zeros(num * num2)
    Iin = torch.tensor([[[1], [0], [0], [0]]], dtype=torch.float64)
    # 调整不同的理论角度
    angle_R1 = torch.tensor([[[np.pi / 4]]])  # np.pi/4 3*np.pi/4
    angle_R2 = torch.tensor([[[0]]])  # 0 np.pi/2
    angle_A = torch.tensor([[[np.pi / 4]]])  # np.pi/4 3*np.pi/4
    start = time.time()

    for epoch in range(num):
        for iter_num in range(num2):
            optimizer = optimizer
            optimizer.zero_grad()
            output = myNet(Y)
            output1 = output[:, :, 0:3 * series_num]
            output2 = output[:, :, 3 * series_num:3 * series_num + 4]
            output3 = output[:, :, 3 * series_num + 4:]

            # 预测可能存在的方位角误差（P-R1-R2-A）
            if align:
                beta1 = output2[:1, :1, 0]
                alpha1 = output2[:1, :1, 1]
                beta2 = output2[:1, :1, 2]
                alpha2 = output2[:1, :1, 3]
            else:
                beta1 = torch.tensor(0)
                alpha1 = torch.tensor(0)
                beta2 = torch.tensor(0)
                alpha2 = torch.tensor(0)

            # 预测波片温度控制系数
            if temperature:
                T = (output3[:1, :1, 0] + 1) * 20 - 20 + 24
                d_delta2 = RetardanceT(wvl, T) * 10
                d_delta1 = 3 * d_delta2
            else:
                d_delta1 = 0
                d_delta2 = 0
                T = torch.tensor(24)

            temp1 = torch.cat([output1[:1, :1, 0:series_num], zero_col], dim=2)
            temp2 = torch.cat([output1[:1, :1, series_num:2 * series_num], zero_col], dim=2)
            temp3 = torch.cat([output1[:1, :1, 2 * series_num:3 * series_num], zero_col], dim=2)

            N_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp1])
            C_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp2])
            S_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp3])

            M_s = Sample(N_g, C_g, S_g)
            M_P = polarizer(beta1)
            M_R1 = waveplate(alpha1 + angle_R1, d1 + d_delta1)
            M_R2 = waveplate(beta2 + angle_R2, d2 + d_delta2)
            M_A = polarizer(alpha2 + angle_A)

            I_g_temp = M_A @ M_R2 @ M_s @ M_R1 @ M_P @ Iin
            I_g_extract = I_g_temp[:, 0, 0]
            I_g = I_g_extract.unsqueeze(0).unsqueeze(0)

            error = (Y[:1, :1, 1:-1] - I_g[:1, :1, 1:-1])

            loss = (torch.norm(error[:1, :1, int(range1 * sample_points):int(range2 * sample_points)], p=2,
                               dim=2).squeeze(0)
                    + beta * torch.norm(output1, p=1, dim=2).squeeze(0)
                    + gama * torch.norm(N_g ** 2 + C_g ** 2 +
                                        S_g ** 2 - 1, p=2, dim=2).squeeze(0))

            loss.backward()
            optimizer.step()
            running_loss[epoch * num2 + iter_num] = loss.item()

    end = time.time()
    time_cost = end - start
    print('iter num: %d, loss=%f' % (epoch * num2 + iter_num + 1, loss.item()), '\n')
    print('time cost: ', time_cost, 's \n')

    return N_g, C_g, S_g, I_g, beta1, alpha1, beta2, alpha2, running_loss, time_cost, T


class CNN(nn.Module):
    def __init__(self, input_length, output_size):
        """
        :param input_length: 输入张量的长度（例如 700）
        :param output_size: 输出张量的目标大小（例如 3 * series_num + 5）
        """
        super(CNN, self).__init__()
        self.conv = nn.Sequential(
            # 第一层卷积：输入通道数为1，输出通道数为64，卷积核大小为3，填充为1
            nn.Conv1d(1, 64, kernel_size=3, padding=1, dtype=torch.float64),
            nn.ReLU(),
            # 第一层池化：将长度减半
            nn.MaxPool1d(kernel_size=2, stride=2),

            # 第二层卷积：输入通道数为64，输出通道数为128，卷积核大小为3，填充为1
            nn.Conv1d(64, 128, kernel_size=3, padding=1, dtype=torch.float64),
            nn.ReLU(),
            # 第二层池化：再次将长度减半
            nn.MaxPool1d(kernel_size=2, stride=2),

            # 展平操作
            nn.Flatten()
        )

        # 计算最终卷积层输出的特征维度
        # 输入长度经过两次池化后变为 input_length // 4
        final_conv_length = input_length // 4
        self.fc = nn.Linear(128 * final_conv_length, output_size, dtype=torch.float64)

    def forward(self, x):
        """
        :param x: 输入张量，形状为 [batch_size, 1, input_length]
        :return: 输出张量，形状为 [1, 1, output_size]
        """
        # 调整输入形状为 (batch_size, channels, length)
        # x = x.permute(0, 2, 1)  # [batch_size, 1, input_length] -> [batch_size, input_length, 1]
        x = self.conv(x)  # 卷积和池化操作
        x = self.fc(x)  # 全连接层
        # 确保输出形状为 [1, 1, output_size]
        return x.view(1, 1, -1)


def SFCNet_CNN(wvl, M_basis, Y, d1, d2, num=2, num2=400,
               f_set=0.05, beta=0.01, gama=1, seed=5, range1=0.05, range2=0.95, align=False, temperature=False):
    # 不分组
    series_num = int(M_basis.shape[2] * f_set)
    zero_col = torch.zeros([1, 1, M_basis.shape[2] - series_num])
    sample_points = int(Y.shape[2])
    setup_seed(seed)
    myNet = CNN(sample_points, 3 * series_num + 5)
    start_lr = 0.001
    params = myNet.parameters()
    optimizer = torch.optim.Adam(params, lr=start_lr)
    running_loss = np.zeros(num * num2)
    Iin = torch.tensor([[[1], [0], [0], [0]]], dtype=torch.float64)
    # 调整不同的理论角度
    angle_R1 = torch.tensor([[[np.pi / 4]]])  # np.pi/4 3*np.pi/4
    angle_R2 = torch.tensor([[[0]]])  # 0 np.pi/2
    angle_A = torch.tensor([[[np.pi / 4]]])  # np.pi/4 3*np.pi/4
    start = time.time()

    for epoch in range(num):
        for iter_num in range(num2):
            optimizer = optimizer
            optimizer.zero_grad()
            output = myNet(Y)
            output1 = output[:, :, 0:3 * series_num]
            output2 = output[:, :, 3 * series_num:3 * series_num + 4]
            output3 = output[:, :, 3 * series_num + 4:]

            # 预测可能存在的方位角误差（P-R1-R2-A）
            if align:
                beta1 = output2[:1, :1, 0]
                alpha1 = output2[:1, :1, 1]
                beta2 = output2[:1, :1, 2]
                alpha2 = output2[:1, :1, 3]
            else:
                beta1 = torch.tensor(0)
                alpha1 = torch.tensor(0)
                beta2 = torch.tensor(0)
                alpha2 = torch.tensor(0)

            # 预测波片温度控制系数
            if temperature:
                T = (output3[:1, :1, 0] + 1) * 20 - 20 + 24
                d_delta2 = RetardanceT(wvl, T) * 10
                d_delta1 = 3 * d_delta2
            else:
                d_delta1 = 0
                d_delta2 = 0
                T = torch.tensor(24)

            temp1 = torch.cat([output1[:1, :1, 0:series_num], zero_col], dim=2)
            temp2 = torch.cat([output1[:1, :1, series_num:2 * series_num], zero_col], dim=2)
            temp3 = torch.cat([output1[:1, :1, 2 * series_num:3 * series_num], zero_col], dim=2)

            N_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp1])
            C_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp2])
            S_g = torch.einsum('ijk,ilk->ilj', [M_basis, temp3])

            M_s = Sample(N_g, C_g, S_g)
            M_P = polarizer(beta1)
            M_R1 = waveplate(alpha1 + angle_R1, d1 + d_delta1)
            M_R2 = waveplate(beta2 + angle_R2, d2 + d_delta2)
            M_A = polarizer(alpha2 + angle_A)

            I_g_temp = M_A @ M_R2 @ M_s @ M_R1 @ M_P @ Iin
            I_g_extract = I_g_temp[:, 0, 0]
            I_g = I_g_extract.unsqueeze(0).unsqueeze(0)

            error = (Y[:1, :1, 1:-1] - I_g[:1, :1, 1:-1])

            loss = (torch.norm(error[:1, :1, int(range1 * sample_points):int(range2 * sample_points)], p=2,
                               dim=2).squeeze(0)
                    + beta * torch.norm(output1, p=1, dim=2).squeeze(0)
                    + gama * torch.norm(N_g ** 2 + C_g ** 2 +
                                        S_g ** 2 - 1, p=2, dim=2).squeeze(0))

            loss.backward()
            optimizer.step()
            running_loss[epoch * num2 + iter_num] = loss.item()
            print('iter num: %d, loss=%f' % (epoch * num2 + iter_num + 1, loss.item()), '\n')

    end = time.time()
    time_cost = end - start
    print('iter num: %d, loss=%f' % (epoch * num2 + iter_num + 1, loss.item()), '\n')
    print('time cost: ', time_cost, 's \n')

    return N_g, C_g, S_g, I_g, beta1, alpha1, beta2, alpha2, running_loss, time_cost, T