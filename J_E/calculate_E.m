function E = calculate_E(Ge,polar)

%% ============================计算电场============================
[~,~,w,m,~,u,~,~,~] = define_E();
% 首先定义电场 幅度为1 理想点源 极化方向为polar
% 由于要计算E的三个分量 所以将球坐标转为笛卡尔坐标
[Jx,Jy,Jz] = sph2cart(polar(1),polar(2),polar(3)); 
J = [Jx,Jy,Jz];
E = zeros(3,1);
E(:,:) = 1i * w * u(m) * Ge * J';


end