function [N,f,w,m,n,u,e,d,k] = define_H()
%% ============================常量定义============================
N = 7; % 介质层数, 应当大于3 否则得修改代码
c = 3e8;
u0 = 4*pi*1e-7; % 磁导率
e0 = 8.8542e-12; % 电导率
f = 300e6; % 频率
w = f * 2 * pi; % 角频率

% 场 源定义，最上层为N层，最底层为1层
m = 2; % 源层数
n = 5; % 场层数

% 每一层的介电常数\磁导率\厚度(分层位置)\波数
u = zeros(N,1);
e = zeros(N,1);
d = zeros(N,1);
k = zeros(N,1);
e(1) = e0;          u(1) = u0;          d(1) = 0;     k(1) = w * sqrt(e(1) * u(1)); % d(1)无意义
e(2) = e0 * 2.6;    u(2) = u0;          d(2) = -1.5;  k(2) = w * sqrt(e(2) * u(2));
e(3) = e0 * 6.5;    u(3) = u0 * 3.2;    d(3) = -1.3;  k(3) = w * sqrt(e(3) * u(3));
e(4) = e0 * 4.2;    u(4) = u0 * 6;      d(4) = -1;    k(4) = w * sqrt(e(4) * u(4));
e(5) = e0 * 6.5;    u(5) = u0 * 3.2;    d(5) = -0.5;  k(5) = w * sqrt(e(5) * u(5));
e(6) = e0 * 2.6;    u(6) = u0;          d(6) = -0.2;  k(6) = w * sqrt(e(6) * u(6));
e(7) = e0;          u(7) = u0;          d(7) = 0;     k(7) = w * sqrt(e(7) * u(7));

temp = u;
u = -e;
e = -temp; %根据对偶原理进行置换

end