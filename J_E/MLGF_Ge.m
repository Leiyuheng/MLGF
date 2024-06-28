function Ge = MLGF_Ge(position_s, position_f)

%% ============================常量定义============================
[N,~,w,m,n,u,e,d,k] = define_E();

x_s = position_s(1); y_s = position_s(2); z_s = position_s(3); % 源位置
x_f = position_f(1); y_f = position_f(2); z_f = position_f(3); % 场位置

% 传递波数定义
k_iz = @(kz, i) sqrt(k(i)^2 - kz.^2); % z代表场位置
k_nm = w * sqrt(e(n) * u(m));

%% ============================计算广义反射系数R============================
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

%% ============================计算透射系数Tmn n>m情况============================
% 首先定义S_j-1,j = S{j}
S_TE = cell(1, N);
S_TM = cell(1, N);
for i = 2:N-1
    S_TM{i} = @(kz) T_TM(kz, i+1, i) ./ (1 - R_TM(kz, i, i-1) .* R_Gen_TM_up{i}(kz) .* exp(1i .* k_iz(kz, i) * 2 * (d(i+1) - d(i))));
    S_TE{i} = @(kz) T_TE(kz, i+1, i) ./ (1 - R_TE(kz, i, i-1) .* R_Gen_TE_up{i}(kz) .* exp(1i .* k_iz(kz, i) * 2 * (d(i+1) - d(i))));
end

% m为源的层，n为场的层
T_Gen_TE_mn = @(kz) 1;
T_Gen_TM_mn = @(kz) 1;
for i = m+1:n-1
    temp_TE = T_Gen_TE_mn;
    func_temp_TE = @(kz) exp(1i .* k_iz(kz, i) * (d(i+1) - d(i))) .* S_TE{i}(kz);
    T_Gen_TE_mn = @(kz) temp_TE(kz) .* func_temp_TE(kz);

    temp_TM = T_Gen_TM_mn;
    func_temp_TM = @(kz) exp(1i .* k_iz(kz, i) * (d(i+1) - d(i))) .* S_TM{i}(kz);
    T_Gen_TM_mn = @(kz) temp_TM(kz) .* func_temp_TM(kz);
end
T_Gen_TM_mn = @(kz) T_Gen_TM_mn(kz) .* S_TM{n}(kz);
T_Gen_TE_mn = @(kz) T_Gen_TE_mn(kz) .* S_TE{n}(kz);

%% ============================计算I_mn============================
% 计算M_m
M_TE_m = @(kz) 1 ./ (1 - R_Gen_TE_down{m}(kz) .* R_Gen_TE_up{m}(kz) .* exp(2 * 1i .* k_iz(kz, m) * (d(m+1) - d(m))));
M_TM_m = @(kz) 1 ./ (1 - R_Gen_TM_down{m}(kz) .* R_Gen_TM_up{m}(kz) .* exp(2 * 1i .* k_iz(kz, m) * (d(m+1) - d(m))));

% 计算I_mn
I_TE_mn = @(kz) M_TE_m(kz) .* T_Gen_TE_mn(kz);
I_TM_mn = @(kz) M_TM_m(kz) .* T_Gen_TM_mn(kz);

%% ============================计算fv与fr及其偏导============================
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

%% ============================计算F(z,z')============================
F_TE = @(kz) fv_TE_zf(kz) .* I_TE_mn(kz) .* fr_TE_zs(kz);
F_TM = @(kz) fv_TM_zf(kz) .* I_TM_mn(kz) .* fr_TM_zs(kz);

F_TE_diff_zfs = @(kz) fv_TE_zf_diff(kz) .* I_TE_mn(kz) .* fr_TE_zs_diff(kz);
F_TM_diff_zfs = @(kz) fv_TM_zf_diff(kz) .* I_TM_mn(kz) .* fr_TM_zs_diff(kz);
F_TM_diff_zf = @(kz) fv_TM_zf_diff(kz) .* I_TM_mn(kz) .* fr_TM_zs(kz);
F_TM_diff_zs = @(kz) fv_TM_zf(kz) .* I_TM_mn(kz) .* fr_TM_zs_diff(kz);
%% ============================计算G============================
G_TE = zeros(3,3,length(x_f));
G_TM = zeros(3,3,length(x_f));
Ge = zeros(3,3,length(x_f));

inf_integral = 300; % 积分路径上限

for i = 1:length(x_f)
    x = x_f(i); y = y_f; z = z_f;
    [phi, rho] = cart2pol(x-x_s, y-y_s);

    f_TE_s0 = @(kz) F_TE(kz) ./ k_iz(kz, m) .* besselj(0, kz * rho) .* kz;  % TE波0阶Sommerfield积分被积函数
    f_TE_s1 = @(kz) F_TE(kz) ./ k_iz(kz, m) .* besselj(1, kz * rho);        % TE波1阶Sommerfield积分被积函数
    
    f_TM_xx_s0 = @(kz) F_TM_diff_zfs(kz) ./ k_iz(kz, m) .* besselj(0, kz * rho) .* kz;      % TMxx波0阶Sommerfield积分被积函数
    f_TM_xx_s1 = @(kz) F_TM_diff_zfs(kz) ./ k_iz(kz, m) .* besselj(1, kz * rho);            % TMxx波1阶Sommerfield积分被积函数
    f_TM_xz_s1 = @(kz) F_TM_diff_zf(kz) ./ k_iz(kz, m) .* besselj(1, kz * rho) .* kz.^2;    % TMxz波1阶Sommerfield积分被积函数
    f_TM_zx_s1 = @(kz) F_TM_diff_zs(kz) ./ k_iz(kz, m) .* besselj(1, kz * rho) .* kz.^2;    % TMxz波1阶Sommerfield积分被积函数
    f_TM_zz_s0 = @(kz) F_TM(kz) ./ k_iz(kz, m) .* besselj(0, kz * rho) .* kz.^3;            % TMxz波1阶Sommerfield积分被积函数
    
    S0_TE = (quadgk(f_TE_s0, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TE_s0, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % 0阶Sommerfield积分
    S1_TE = (quadgk(f_TE_s1, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TE_s1, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % 1阶Sommerfield积分
    
    S0_TM_xx = (quadgk(f_TM_xx_s0, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TM_xx_s0, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % TM波xx分量0阶Sommerfield积分
    S1_TM_xx = (quadgk(f_TM_xx_s1, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TM_xx_s1, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % TM波xx分量1阶Sommerfield积分
    S1_TM_xz = (quadgk(f_TM_xz_s1, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TM_xz_s1, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % TM波xx分量0阶Sommerfield积分
    S1_TM_zx = (quadgk(f_TM_zx_s1, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TM_zx_s1, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % TM波xx分量1阶Sommerfield积分
    S0_TM_zz = (quadgk(f_TM_zz_s0, 0, -0.1i,'MaxIntervalCount',1e5)+quadgk(f_TM_zz_s0, -0.1i, inf_integral-0.1i,'MaxIntervalCount',1e5)) /(2*pi);  % TM波xx分量0阶Sommerfield积分

    G_TE_xx = 1i / (2*rho) * cos(2*phi) * S1_TE + 1i / 4 * (1-cos(2*phi)) * S0_TE;
    G_TE_xy = 1i / (2*rho) * sin(2*phi) * S1_TE - 1i / 4 * sin(2*phi) * S0_TE;
    G_TE_yx = G_TE_xy;
    G_TE_yy = -1i / (2*rho) * cos(2*phi) * S1_TE + 1i / 4 * (1+cos(2*phi)) * S0_TE;

    G_TM_xx = -1i / (2*rho) * cos(2*phi) * S1_TM_xx + 1i / 4 * (1+cos(2*phi)) * S0_TM_xx;
    G_TM_xy = -1i / (2*rho) * sin(2*phi) * S1_TM_xx + 1i / 4 * sin(2*phi) * S0_TM_xx;
    G_TM_xz = -1i / 2 * cos(phi) * S1_TM_xz;
    G_TM_yx = G_TM_xy;
    G_TM_yy = 1i / (2*rho) * cos(2*phi) * S1_TM_xx + 1i / 4 * (1-cos(2*phi)) * S0_TM_xx;
    G_TM_yz = -1i / 2 * sin(phi) * S1_TM_xz;
    G_TM_zx = 1i / 2 * cos(phi) * S1_TM_zx;
    G_TM_zy = 1i / 2 * sin(phi) * S1_TM_zx;
    G_TM_zz = 1i / 2 * S0_TM_zz;


    G_TE(:,:,i) = [G_TE_xx G_TE_xy 0; G_TE_yx G_TE_yy 0; 0 0 0];
    G_TM(:,:,i) = [G_TM_xx G_TM_xy G_TM_xz; G_TM_yx G_TM_yy G_TM_yz; G_TM_zx G_TM_zy G_TM_zz];
    Ge(:,:,i) = G_TE(:,:,i) + G_TM(:,:,i)/k_nm^2;
end

% 保存数据
save('G_TE.mat','G_TE');
save('G_TM.mat','G_TM');
save('Ge.mat','Ge');


end





