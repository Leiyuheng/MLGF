function H = calculate_H(Gh,polar)

%% ============================计算磁场============================
[~,~,w,m,~,u,e,~,~] = define_H();
% 首先定义电场 幅度为1 理想点源 极化方向为polar
% 由于要计算E的三个分量 所以将球坐标转为笛卡尔坐标
[Jx,Jy,Jz] = sph2cart(polar(1),polar(2),polar(3));
J = [Jx,Jy,Jz];
M = J;
H = zeros(3,1);
H(:,:) = 1i * w * (-u(m)) * Gh * M'; % define中由于对偶原理更换了e和u 这里要换回来


end