close all
clear
clc

%% ======================常量定义=======================
N = 7; % 介质层数, 应当大于3 否则得修改代码
u0 = 4 * pi * 1e-7; % 磁导率
e0 = 8.8542e-12; % 电导率
f = 300e6; % 频率
w = f * 2 * pi; % 角频率

% 场 源定义，最上层为N层，最底层为1层
m = 2; % 源层数
n = 5; % 场层数
x_s = 0; y_s = 0; z_s = -1.4; % 源位置
x_f = -3:0.1:3; y_f = 1; z_f = -0.3; % 场位置

% 每一层的介电常数\磁导率\厚度(分层位置)\波数
u = zeros(N,1);
e = zeros(N,1);
d = zeros(N,1);
k = zeros(N,1);
e(1) = e0;    u(1) = u0;    d(1) = 0;     k(1) = w * sqrt(e(1) * u(1)); % d(0)无意义
e(2) = e0 * 2.6; u(2) = u0;    d(2) = -1.5;  k(2) = w * sqrt(e(2) * u(2));
e(3) = e0 * 6.5; u(3) = u0 * 3.2; d(3) = -1.3;  k(3) = w * sqrt(e(3) * u(3));
e(4) = e0 * 4.2; u(4) = u0 * 6;  d(4) = -1;    k(4) = w * sqrt(e(4) * u(4));
e(5) = e0 * 6.5; u(5) = u0 * 3.2; d(5) = -0.5;  k(5) = w * sqrt(e(5) * u(5));
e(6) = e0 * 2.6; u(6) = u0;    d(6) = -0.2;  k(6) = w * sqrt(e(6) * u(6));
e(7) = e0;    u(7) = u0;    d(7) = 0;     k(7) = w * sqrt(e(7) * u(7));
% 传递波数定义
k_iz = @(kz, i) sqrt(k(i)^2 - kz.^2); % z代表场位置

%% ===============================计算广义反射系数R============================
% 首先定义反射系数与透射系数
% 广义反射系数可以一层一层的递推出R_i,i+1
R_TM = @(kz, p, q) (e(q) * k_iz(kz, p) - e(p) * k_iz(kz, q)) ./ (e(q) * k_iz(kz, p) + e(p) * k_iz(kz, q)); % p入射到q的反射系数
R_TE = @(kz, p, q) (u(q) * k_iz(kz, p) - u(p) * k_iz(kz, q)) ./ (u(q) * k_iz(kz, p) + u(p) * k_iz(kz, q));
T_TM = @(kz, p, q) R_TM(kz, p, q) + 1; % p入射到q的透射系数
T_TE = @(kz, p, q) R_TE(kz, p, q) + 1;

% 边界条件 最后一层广义反射系数R_n,n+1等于普通反射系数

% 向上的广义反射系数,N层N-1个边界
R_Gen_TM_up = cell(1, N); % 1层到6层有意义
R_Gen_TE_up = cell(1, N);

R_Gen_TM_up{N} = @(kz) 0;
R_Gen_TE_up{N} = @(kz) 0;

R_Gen_TM_up{N-1} = @(kz) R_TM(kz, N-1, N); % GEN_R_TM(i,i+1)
R_Gen_TE_up{N-1} = @(kz) R_TE(kz, N-1, N);

for i = N-2:-1:1
    R_Gen_TM_up{i} = @(kz) (R_TM(kz, i, i+1) + R_Gen_TM_up{i+1}(kz) .* exp(2 * 1i .* k_iz(kz, i+1) * (d(i+2) - d(i+1)))) ./ ...
                         (1 + R_TM(kz, i, i+1) .* R_Gen_TM_up{i+1}(kz) .* exp(2 * 1i .* k_iz(kz, i+1) * (d(i+2) - d(i+1))));
    R_Gen_TE_up{i} = @(kz) (R_TE(kz, i, i+1) + R_Gen_TE_up{i+1}(kz) .* exp(2 * 1i .* k_iz(kz, i+1) * (d(i+2) - d(i+1)))) ./ ...
                         (1 + R_TE(kz, i, i+1) .* R_Gen_TE_up{i+1}(kz) .* exp(2 * 1i .* k_iz(kz, i+1) * (d(i+2) - d(i+1)))); % 注意d在这里应该是i+1的层厚度，向上入射
end

% 向下的广义反射系数
R_Gen_TM_down = cell(1, N); % 2层到7层有意义
R_Gen_TE_down = cell(1, N);

R_Gen_TM_down{1} = @(kz) 0;
R_Gen_TE_down{1} = @(kz) 0;

R_Gen_TM_down{2} = @(kz) R_TM(kz, 2, 1); % GEN_R_TM(i,i-1)
R_Gen_TE_down{2} = @(kz) R_TE(kz, 2, 1);

for i = 3:N
    R_Gen_TM_down{i} = @(kz) (R_TM(kz, i, i-1) + R_Gen_TM_down{i-1}(kz) .* exp(2 * 1i .* k_iz(kz, i-1) * (d(i) - d(i-1)))) ./ ...
                         (1 + R_TM(kz, i, i-1) .* R_Gen_TM_down{i-1}(kz) .* exp(2 * 1i .* k_iz(kz, i-1) * (d(i) - d(i-1))));
    R_Gen_TE_down{i} = @(kz) (R_TE(kz, i, i-1) + R_Gen_TE_down{i-1}(kz) .* exp(2 * 1i .* k_iz(kz, i-1) * (d(i) - d(i-1)))) ./ ...
                         (1 + R_TE(kz, i, i-1) .* R_Gen_TE_down{i-1}(kz) .* exp(2 * 1i .* k_iz(kz, i-1) * (d(i) - d(i-1)))); % 注意d在这里应该是i-1的层厚度，向下入射
end

%% ================================计算透射系数Tmn n>m情况==============================
% 首先定义S_j-1,j = S{j}
S_TE = cell(1, N);
S_TM = cell(1, N);
for j = 2:N-1
    S_TM{j} = @(kz) T_TM(kz, j+1, j) ./ (1 - R_TM(kz, j, j-1) .* R_Gen_TM_up{j}(kz) .* exp(1i .* k_iz(kz, j) * 2 * (d(j+1) - d(j))));
    S_TE{j} = @(kz) T_TE(kz, j+1, j) ./ (1 - R_TE(kz, j, j-1) .* R_Gen_TE_up{j}(kz) .* exp(1i .* k_iz(kz, j) * 2 * (d(j+1) - d(j))));
end

% m为源的层，n为场的层
T_Gen_TE_mn = @(kz) 1;
T_Gen_TM_mn = @(kz) 1;
for j = m+1:n-1
    temp_TE = T_Gen_TE_mn;
    func_temp_TE = @(kz) exp(1i .* k_iz(kz, j) * (d(j+1) - d(j))) .* S_TE{j}(kz);
    T_Gen_TE_mn = @(kz) temp_TE(kz) .* func_temp_TE(kz);

    temp_TM = T_Gen_TM_mn;
    func_temp_TM = @(kz) exp(1i .* k_iz(kz, j) * (d(j+1) - d(j))) .* S_TM{j}(kz);
    T_Gen_TM_mn = @(kz) temp_TM(kz) .* func_temp_TM(kz);
end
T_Gen_TM_mn = @(kz) T_Gen_TM_mn(kz) .* S_TM{n}(kz);
T_Gen_TE_mn = @(kz) T_Gen_TE_mn(kz) .* S_TE{n}(kz);

%% ======================================计算I_mn======================================
% 计算M_m
M_TE_m = @(kz) 1 ./ (1 - R_Gen_TE_down{m}(kz) .* R_Gen_TE_up{m}(kz) .* exp(2 * 1i .* k_iz(kz, m) * (d(m+1) - d(m))));
M_TM_m = @(kz) 1 ./ (1 - R_Gen_TM_down{m}(kz) .* R_Gen_TM_up{m}(kz) .* exp(2 * 1i .* k_iz(kz, m) * (d(m+1) - d(m))));

% 计算I_mn
I_TE_mn = @(kz) M_TE_m(kz) .* T_Gen_TE_mn(kz);
I_TM_mn = @(kz) M_TM_m(kz) .* T_Gen_TM_mn(kz);

%% ====================================计算fv与fr及其偏导======================================
% z_f(即z)是场的位置 z_s(即z')是源的位置
fv_TM_zf = @(kz) exp(1i .* k_iz(kz, n) .* (z_f - d(n))) + R_Gen_TM_up{n}(kz) .* exp(1i .* k_iz(kz, n) .* (2 * d(n+1) - d(n) - z_f));
fv_TE_zf = @(kz) exp(1i .* k_iz(kz, n) .* (z_f - d(n))) + R_Gen_TE_up{n}(kz) .* exp(1i .* k_iz(kz, n) .* (2 * d(n+1) - d(n) - z_f));

fr_TM_zs = @(kz) exp(1i .* k_iz(kz, m) * (d(m+1) - z_s)) + R_Gen_TM_down{m}(kz) .* exp(1i .* k_iz(kz, m) .* (d(m+1) - 2 * d(m) + z_s));
fr_TE_zs = @(kz) exp(1i .* k_iz(kz, m) * (d(m+1) - z_s)) + R_Gen_TE_down{m}(kz) .* exp(1i .* k_iz(kz, m) .* (d(m+1) - 2 * d(m) + z_s));

% 计算fv与fr偏导数
fv_TM_zf_diff = @(kz) 1i .* k_iz(kz, n) .* (exp(1i .* k_iz(kz, n) .* (z_f - d(n))) - R_Gen_TM_up{n}(kz) .* exp(1i .* k_iz(kz, n) .* (2 * d(n+1) - d(n) - z_f)));
fv_TE_zf_diff = @(kz) 1i .* k_iz(kz, n) .* (exp(1i .* k_iz(kz, n) * (z_f - d(n))) - R_Gen_TE_up{n}(kz) .* exp(1i .* k_iz(kz, n) .* (2 * d(n+1) - d(n) - z_f)));

fr_TM_zs_diff = @(kz) -1i .* k_iz(kz, m) .* (exp(1i .* k_iz(kz, m) .* (d(m+1) - z_s)) - R_Gen_TM_down{m}(kz) .* exp(1i .* k_iz(kz, m) .* (d(m+1) - 2 * d(m) + z_s)));
fr_TE_zs_diff = @(kz) -1i .* k_iz(kz, m) .* (exp(1i .* k_iz(kz, m) .* (d(m+1) - z_s)) - R_Gen_TE_down{m}(kz) .* exp(1i .* k_iz(kz, m) .* (d(m+1) - 2 * d(m) + z_s)));

%% ======================================计算F(z,z')=============================================
F_TE = @(kz) fv_TE_zf(kz) .* I_TE_mn(kz) .* fr_TE_zs(kz);
F_TM = @(kz) fv_TM_zf(kz) .* I_TM_mn(kz) .* fr_TM_zs(kz);

%% =====================================计算GF===================================================
% 首先计算TE的GF
x = 0.1;
y = 1;
rho = sqrt(x^2 + y^2);
tilde_f = @(kz) F_TE(kz) ./ k_iz(kz, m); % 被积函数
integrand = @(kz) (tilde_f(kz) .* besselj(0, kz * rho) .* kz);
S0_TE = (1 / (2 * pi)) * quadgk(integrand, 0, Inf);


