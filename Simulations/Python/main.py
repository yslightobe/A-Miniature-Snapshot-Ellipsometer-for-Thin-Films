from NetParts import *

print('main 成功运行 \n')
# 读取测量数据
current_dir = os.path.dirname(__file__)
root = os.path.abspath(os.path.join(current_dir, '..'))
file_path = root + r"\Matlab\data.csv"
d1_np = read_col(file_path, 'd1')
d2_np = read_col(file_path, 'd2')
wvl_um_np = read_col(file_path, 'wvl_um')
hyper = torch.from_numpy(np.array(pd.read_csv(root + r"\Matlab\hyper.csv"), dtype=np.float64))
Y_np = read_col(file_path, 'Y')
M_basis_np = np.array(pd.read_csv
                      (root + r"\Matlab\M_basis.csv"),
                      dtype=np.float64)
# 将数据转换为GPU张量
d1 = torch.from_numpy(d1_np).unsqueeze(0)
d2 = torch.from_numpy(d2_np).unsqueeze(0)
d1 = d1 % 360 * np.pi / 180
d2 = d2 % 360 * np.pi / 180
wvl_um = torch.from_numpy(wvl_um_np).unsqueeze(0)
Y = torch.from_numpy(Y_np).unsqueeze(0)
M_basis = torch.from_numpy(M_basis_np).unsqueeze(0)

N_g, C_g, S_g, I_g, beta1, alpha1, beta2, alpha2, running_loss, time_cost, T = SFCNet(wvl_um, M_basis, Y, d1, d2, num=2,
                                                                              num2=100, f_set=hyper[0, 0],
                                                                              beta=hyper[0, 1],
                                                                              gama=hyper[0, 2], seed=5, range1=0,
                                                                              range2=1, align=bool(hyper[0, 3]), temperature=bool(hyper[0, 4]))

# 储存数据
out_N = N_g.squeeze(0).squeeze(0).detach().numpy()
out_C = C_g.squeeze(0).squeeze(0).detach().numpy()
out_S = S_g.squeeze(0).squeeze(0).detach().numpy()
out_I = I_g.squeeze(0).squeeze(0).detach().numpy()
out_b1 = beta1.squeeze(0).squeeze(0).detach().numpy() * 180 / np.pi
out_a1 = alpha1.squeeze(0).squeeze(0).detach().numpy() * 180 / np.pi
out_b2 = beta2.squeeze(0).squeeze(0).detach().numpy() * 180 / np.pi
out_a2 = alpha2.squeeze(0).squeeze(0).detach().numpy() * 180 / np.pi
out_beta2 = beta2.squeeze(0).squeeze(0).detach().numpy() * 180 / np.pi
x = pd.DataFrame({'out_N': out_N, 'out_C': out_C,
                  'out_S': out_S, 'I_predict': out_I})
x.to_csv(root + r"\Matlab\result.csv")
running_loss = np.append(running_loss, time_cost)
running_loss = np.append(running_loss, out_b1)
running_loss = np.append(running_loss, out_a1)
running_loss = np.append(running_loss, out_b2)
running_loss = np.append(running_loss, out_a2)
running_loss = np.append(running_loss, T.squeeze(0).squeeze(0).detach().numpy())
y = pd.DataFrame({'loss_time_eps': running_loss})
y.to_csv(root + r"\Matlab\loss.csv")
