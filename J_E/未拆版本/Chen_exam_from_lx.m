clear
N=7; %介质层数
j=1j;
u0=4*pi*1e-7;
e0=8.8542e-12;
f=300e6;
w=f*2*pi;
%%==============介电常数定义==============
e=zeros(N,1);
e(1)=e0;
e(2)=e0*2.6;
e(3)=e0*6.5;
e(4)=e0*4.2;
e(5)=e0*6.5;
e(6)=e0*2.6;
e(7)=e0;
%%==================磁导率定义=============
u=zeros(N,1);
u(1)=u0;
u(2)=u0;
u(3)=u0*3.2;
u(4)=u0*6;
u(5)=u0*3.2;
u(6)=u0;
u(7)=u0;

G_TE=zeros(3,3,61);
G_TM=zeros(3,3,61);
J=zeros(3,1);
k=zeros(N,1);
for i=1:N
    k(i)=w*sqrt(e(i)*u(i));  % 每一层的波数
end
d=zeros(N,1);
d(2)=-1.5;d(3)=-1.3;d(4)=-1;d(5)=-0.5;d(6)=-0.2;d(7)=0;
m=2;n=5;  %%场在n层，源在m层
z_f=-0.3;  %%没有将z_f设为自变量 S观察点位置
z_s=-1.4;x_s=0;y_s=0;  %%源点位置
ro=1;   %%定义ρ
phi=@(x_f,y_f) acos((x_f-x_s)./sqrt((x_f-x_s).^2+(y_f-y_s).^2));
kz=@(kro,i) sqrt(k(i)^2-kro.^2);
kmn=w*sqrt(e(n)*u(m));

%%======================================定义反射系数和透射系数======================================
R_TM=@(kro,p,q) (e(q).*kz(kro,p)-e(p).*kz(kro,q))./(e(q).*kz(kro,p)+e(p).*kz(kro,q)); % p入射到q的反射系数
R_TE=@(kro,p,q) (u(q).*kz(kro,p)-u(p).*kz(kro,q))./(u(q).*kz(kro,p)+u(p).*kz(kro,q));
T_TM=@(kro,p,q) R_TM(kro,p,q)+1; % p入射到q的透射系数
T_TE=@(kro,p,q) R_TE(kro,p,q)+1;
%%====================================向上的广义反射系数==========================================
gma_r_TM6=@(kro) R_TM(kro,6,7);
gma_r_TM5=@(kro) (R_TM(kro,5,6)+gma_r_TM6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))))./(1+R_TM(kro,5,6).*gma_r_TM6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))));
gma_r_TM4=@(kro) (R_TM(kro,4,5)+gma_r_TM5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))))./(1+R_TM(kro,4,5).*gma_r_TM5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))));
gma_r_TM3=@(kro) (R_TM(kro,3,4)+gma_r_TM4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))))./(1+R_TM(kro,3,4).*gma_r_TM4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))));
gma_r_TM2=@(kro) (R_TM(kro,2,3)+gma_r_TM3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))))./(1+R_TM(kro,2,3).*gma_r_TM3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))));
gma_r_TM1=@(kro) (R_TM(kro,1,2)+gma_r_TM2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))))./(1+R_TM(kro,1,2).*gma_r_TM2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))));

gma_r_TE6=@(kro) R_TE(kro,6,7);
gma_r_TE5=@(kro) (R_TE(kro,5,6)+gma_r_TE6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))))./(1+R_TE(kro,5,6).*gma_r_TE6(kro).*exp(j*kz(kro,6)*2*(d(7)-d(6))));
gma_r_TE4=@(kro) (R_TE(kro,4,5)+gma_r_TE5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))))./(1+R_TE(kro,4,5).*gma_r_TE5(kro).*exp(j*kz(kro,5)*2*(d(6)-d(5))));
gma_r_TE3=@(kro) (R_TE(kro,3,4)+gma_r_TE4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))))./(1+R_TE(kro,3,4).*gma_r_TE4(kro).*exp(j*kz(kro,4)*2*(d(5)-d(4))));
gma_r_TE2=@(kro) (R_TE(kro,2,3)+gma_r_TE3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))))./(1+R_TE(kro,2,3).*gma_r_TE3(kro).*exp(j*kz(kro,3)*2*(d(4)-d(3))));
gma_r_TE1=@(kro) (R_TE(kro,1,2)+gma_r_TE2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))))./(1+R_TE(kro,1,2).*gma_r_TE2(kro).*exp(j*kz(kro,2)*2*(d(3)-d(2))));
%%====================================向下的广义反射系数=========================================
gma_l_TM2=@(kro) R_TM(kro,2,1);
gma_l_TM3=@(kro) (R_TM(kro,3,2)+gma_l_TM2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))))./(1+R_TM(kro,3,2).*gma_l_TM2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))));
gma_l_TM4=@(kro) (R_TM(kro,4,3)+gma_l_TM3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))))./(1+R_TM(kro,4,3).*gma_l_TM3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))));
gma_l_TM5=@(kro) (R_TM(kro,5,4)+gma_l_TM4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))))./(1+R_TM(kro,5,4).*gma_l_TM4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))));
gma_l_TM6=@(kro) (R_TM(kro,6,5)+gma_l_TM5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))))./(1+R_TM(kro,6,5).*gma_l_TM5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))));
gma_l_TM7=@(kro) (R_TM(kro,7,6)+gma_l_TM6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))))./(1+R_TM(kro,7,6).*gma_l_TM6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))));

gma_l_TE2=@(kro) R_TE(kro,2,1);
gma_l_TE3=@(kro) (R_TE(kro,3,2)+gma_l_TE2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))))./(1+R_TE(kro,3,2).*gma_l_TE2(kro).*exp(j.*kz(kro,2)*2*(d(3)-d(2))));
gma_l_TE4=@(kro) (R_TE(kro,4,3)+gma_l_TE3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))))./(1+R_TE(kro,4,3).*gma_l_TE3(kro).*exp(j.*kz(kro,3)*2*(d(4)-d(3))));
gma_l_TE5=@(kro) (R_TE(kro,5,4)+gma_l_TE4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))))./(1+R_TE(kro,5,4).*gma_l_TE4(kro).*exp(j.*kz(kro,4)*2*(d(5)-d(4))));
gma_l_TE6=@(kro) (R_TE(kro,6,5)+gma_l_TE5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))))./(1+R_TE(kro,6,5).*gma_l_TE5(kro).*exp(j.*kz(kro,5)*2*(d(6)-d(5))));
gma_l_TE7=@(kro) (R_TE(kro,7,6)+gma_l_TE6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))))./(1+R_TE(kro,7,6).*gma_l_TE6(kro).*exp(j.*kz(kro,6)*2*(d(7)-d(6))));
%%定义M
M_TM_m=@(kro) (1-gma_l_TM2(kro).*gma_r_TM2(kro).*exp(2*j.*kz(kro,m)*(d(m+1)-d(m)))).^(-1);
M_TE_m=@(kro) (1-gma_l_TE2(kro).*gma_r_TE2(kro).*exp(2*j.*kz(kro,m)*(d(m+1)-d(m)))).^(-1);
%%定义S(kro)
S3_TE=@(kro) T_TE(kro,2,3)./(1-R_TE(kro,3,2).*gma_r_TE3(kro).*exp(j*2.*kz(kro,3)*(d(4)-d(3))));
S4_TE=@(kro) T_TE(kro,3,4)./(1-R_TE(kro,4,3).*gma_r_TE4(kro).*exp(j*2.*kz(kro,4)*(d(5)-d(4))));
S5_TE=@(kro) T_TE(kro,4,5)./(1-R_TE(kro,5,4).*gma_r_TE5(kro).*exp(j*2.*kz(kro,5)*(d(6)-d(5))));

S3_TM=@(kro) T_TM(kro,2,3)./(1-R_TM(kro,3,2).*gma_r_TM3(kro).*exp(j*2.*kz(kro,3)*(d(4)-d(3))));
S4_TM=@(kro) T_TM(kro,3,4)./(1-R_TM(kro,4,3).*gma_r_TM4(kro).*exp(j*2.*kz(kro,4)*(d(5)-d(4))));
S5_TM=@(kro) T_TM(kro,4,5)./(1-R_TM(kro,5,4).*gma_r_TM5(kro).*exp(j*2.*kz(kro,5)*(d(6)-d(5))));
%%定义T25(kro)
T25_TE=@(kro) (exp(j.*kz(kro,3)*(d(4)-d(3))).*S3_TE(kro)).*(exp(j.*kz(kro,4)*(d(5)-d(4))).*S4_TE(kro)).*S5_TE(kro);
T25_TM=@(kro) (exp(j.*kz(kro,3)*(d(4)-d(3))).*S3_TM(kro)).*(exp(j.*kz(kro,4)*(d(5)-d(4))).*S4_TM(kro)).*S5_TM(kro);
%%定义I25(kro)
I25_TE=@(kro) M_TE_m(kro).*T25_TE(kro);
I25_TM=@(kro) M_TM_m(kro).*T25_TM(kro);
%%定义fv(kro)，fr(kro)
fv_TM=@(kro) exp(j.*kz(kro,n)*(z_f-d(n)))+gma_r_TM5(kro).*exp(j.*kz(kro,n)*(2*d(n+1)-d(n)-z_f));
fv_TE=@(kro) exp(j.*kz(kro,n)*(z_f-d(n)))+gma_r_TE5(kro).*exp(j.*kz(kro,n)*(2*d(n+1)-d(n)-z_f));

fr_TM=@(kro) exp(j.*kz(kro,m)*(d(m+1)-z_s))+gma_l_TM2(kro).*exp(j.*kz(kro,m)*(d(m+1)-2*d(m)+z_s));
fr_TE=@(kro) exp(j.*kz(kro,m)*(d(m+1)-z_s))+gma_l_TE2(kro).*exp(j.*kz(kro,m)*(d(m+1)-2*d(m)+z_s));  
%%定义F(kro)
F_TM=@(kro) fr_TM(kro).*fv_TM(kro).*I25_TM(kro);
F_TE=@(kro) fr_TE(kro).*fv_TE(kro).*I25_TE(kro);
% S0xxTE_c=@(kro) F_TE(kro)./sqrt(k(m)^2-kro.^2);
% S0xxTE=quadgk(S0xxTE_c,0+0.001j,260+0.001j)+quadgk(S0xxTE_c,0,0.001j);     %%未进行变量替换

fv_TM_zf=@(kro) j.*kz(kro,n).*(exp(j.*kz(kro,n)*(z_f-d(n)))-gma_r_TM5(kro).*exp(j.*kz(kro,n)*(2*d(n+1)-d(n)-z_f)));  %%fv and fr的导数
fv_TE_zf=@(kro) j.*kz(kro,n).*(exp(j.*kz(kro,n)*(z_f-d(n)))-gma_r_TE5(kro).*exp(j.*kz(kro,n)*(2*d(n+1)-d(n)-z_f)));

fr_TM_zs=@(kro) -j.*kz(kro,m).*(exp(j.*kz(kro,m)*(d(m+1)-z_s))-gma_l_TM2(kro).*exp(j.*kz(kro,m)*(d(m+1)-2*d(m)+z_s)));
fr_TE_zs=@(kro) -j.*kz(kro,m).*(exp(j.*kz(kro,m)*(d(m+1)-z_s))-gma_l_TE2(kro).*exp(j.*kz(kro,m)*(d(m+1)-2*d(m)+z_s)));

%% 计算G_TE部分
y=1;i=1;
for x=-3:0.1:3
    ro=sqrt(x^2+y^2);
    S0xxTE_c=@(krox) F_TE(krox-0.001j)./kz(krox-0.001j,m).*(krox-0.001j).*besselj(0,(krox-0.001j).*ro);  %%Gxx中0阶索末菲积分
    S0xxTE_c2=@(kroy) F_TE(kroy*j)./kz(kroy*j,m).*(kroy*j).*besselj(0,kroy*j.*ro);
    S0xxTE=(quadgk(S0xxTE_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S0xxTE_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;    %%进行了变量替换
    
    S1xxTE_c=@(krox) F_TE(krox-0.001j)./kz(krox-0.001j,m).*besselj(1,(krox-0.001j).*ro);   %%Gxx 1阶索末菲积分
    S1xxTE_c2=@(kroy) F_TE(kroy*j)./kz(kroy*j,m).*besselj(1,kroy*j.*ro);
    S1xxTE=(quadgk(S1xxTE_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S1xxTE_c2,0,-0.001,'MaxIntervalCount',1e5))/(2*pi);
    
    Gxx_TE=@(x_f,y_f) j/2/ro.*cos(2.*phi(x_f,y_f))*S1xxTE+j/4.*(1-cos(2.*phi(x_f,y_f)))*S0xxTE;  %%GeTE矩阵的内容
    Gxy_TE=@(x_f,y_f) j/2/ro.*sin(2.*phi(x_f,y_f))*S1xxTE-j/4.*sin(2.*phi(x_f,y_f))*S0xxTE;
    Gyx_TE=@(x_f,y_f) Gxy_TE(x_f,y_f);
    Gyy_TE=@(x_f,y_f) -j/2/ro.*cos(2*phi(x_f,y_f))*S1xxTE+j/4.*(1+cos(2.*phi(x_f,y_f)))*S0xxTE;
    
    %% 开始计算G_TM部分
    
    S1xxTM_c=@(krox) fv_TM_zf(krox-0.001j).*fr_TM_zs(krox-0.001j).*I25_TM(krox-0.001j)./kz(krox-0.001j,m).*besselj(1,(krox-0.001j).*ro);
    S1xxTM_c2=@(kroy) fv_TM_zf(j*kroy).*fr_TM_zs(j*kroy).*I25_TM(j*kroy)./kz(j*kroy,m).*besselj(1,j*kroy.*ro);
    S1xxTM=(quadgk(S1xxTM_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S1xxTM_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;  
    
    S0xxTM_c=@(krox) fv_TM_zf(krox-0.001j).*fr_TM_zs(krox-0.001j).*I25_TM(krox-0.001j)./kz(krox-0.001j,m).*(krox-0.001j).*besselj(0,(krox-0.001j).*ro);  %%Gxx中0阶索末菲积分
    S0xxTM_c2=@(kroy) fv_TM_zf(j*kroy).*fr_TM_zs(j*kroy).*I25_TM(j*kroy)./kz(kroy*j,m).*(kroy*j).*besselj(0,kroy*j.*ro);
    S0xxTM=(quadgk(S0xxTM_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S0xxTM_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;    
    
    S1xzTM_c=@(krox) fv_TM_zf(krox-0.001j).*fr_TM(krox-0.001j).*I25_TM(krox-0.001j)./kz(krox-0.001j,m).*(krox-0.001j).^2.*besselj(1,(krox-0.001j).*ro);
    S1xzTM_c2=@(kroy) fv_TM_zf(j*kroy).*fr_TM(j*kroy).*I25_TM(j*kroy)./kz(j*kroy,m).*(j*kroy).^2.*besselj(1,(j*kroy).*ro);
    S1xzTM=(quadgk(S1xzTM_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S1xzTM_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;
    
    S1zxTM_c=@(krox) fr_TM_zs(krox-0.001j).*fv_TM(krox-0.001j).*I25_TM(krox-0.001j)./kz(krox-0.001j,m).*(krox-0.001j).^2.*besselj(1,(krox-0.001j).*ro);
    S1zxTM_c2=@(kroy) fr_TM_zs(j*kroy).*fv_TM(j*kroy).*I25_TM(j*kroy)./kz(j*kroy,m).*(j*kroy).^2.*besselj(1,(j*kroy).*ro);
    S1zxTM=(quadgk(S1zxTM_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S1zxTM_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;
    
    S0zzTM_c=@(krox) F_TM(krox-0.001j)./kz(krox-0.001j,m).*(krox-0.001j).^3.*besselj(0,(krox-0.001j).*ro);
    S0zzTM_c2=@(kroy) F_TM(j*kroy)./kz(j*kroy,m).*(j*kroy).^3.*besselj(0,(j*kroy).*ro);
    S0zzTM=(quadgk(S0zzTM_c,0,260,'MaxIntervalCount',1e5)+j*quadgk(S0zzTM_c2,0,-0.001,'MaxIntervalCount',1e5))/2/pi;
    
    Gxx_TM=@(x_f,y_f) -j/ro/2.*cos(2*phi(x_f,y_f))*S1xxTM+j/4.*(1+cos(2.*phi(x_f,y_f)))*S0xxTM;
    Gxy_TM=@(x_f,y_f) -j/ro/2.*sin(2*phi(x_f,y_f))*S1xxTM+j/4.*sin(2*phi(x_f,y_f))*S0xxTM;
    Gxz_TM=@(x_f,y_f) -j/2.*cos(phi(x_f,y_f))*S1xzTM;
    Gyx_TM=@(x_f,y_f) Gxy_TM(x_f,y_f);
    Gyy_TM=@(x_f,y_f) j/ro/2.*cos(2*phi(x_f,y_f))*S1xxTM+j/4.*(1-cos(2*phi(x_f,y_f)))*S0xxTM;
    Gyz_TM=@(x_f,y_f) -j/2.*sin(phi(x_f,y_f))*S1xzTM;
    Gzx_TM=@(x_f,y_f) j/2.*cos(phi(x_f,y_f))*S1zxTM;
    Gzy_TM=@(x_f,y_f) j/2.*sin(phi(x_f,y_f))*S1zxTM;
    Gzz_TM=@(x_f,y_f) j/2.*S0zzTM;
    
    %% 开始计算各点Ge
    
    G_TE(:,:,i)=[Gxx_TE(x,y),Gxy_TE(x,y),0;Gyx_TE(x,y),Gyy_TE(x,y),0;0,0,0];
    G_TM(:,:,i)=[Gxx_TM(x,y),Gxy_TM(x,y),Gxz_TM(x,y);Gyx_TM(x,y),Gyy_TM(x,y),Gyz_TM(x,y);Gzx_TM(x,y),Gzy_TM(x,y),Gzz_TM(x,y)];
    Ge(:,:,i)=G_TE(:,:,i)+G_TM(:,:,i)/(kmn^2);
    i=i+1;
end
save('Ge_liu.mat','Ge');
save('G_TE_liu.mat','G_TE');
save('G_TM_liu.mat','G_TM');
%% 计算电场
[x_d,y_d,z_d]=sph2cart(pi/6,7*pi/18,1);
J=[x_d;y_d;z_d];
E=zeros(3,1,61); 
for i=1:61
    E(:,:,i)=j*w*u(m)*Ge(:,:,i)*J;
end
Ex=zeros(61,1);
for i=1:61
    Ex(i)=abs(E(1,1,i));
    Ex1(i)=E(1,1,i);
end
Ey=zeros(61,1);
for i=1:61
    Ey(i)=abs(E(2,1,i));
    Ey1(i)=E(2,1,i);
end
Ez=zeros(61,1);
for i=1:61
    Ez(i)=abs(E(3,1,i));
    Ez1(i)=E(3,1,i);
end
for i=1:61
    E_total(i)=sqrt((Ex1(i))^2+(Ey1(i))^2+(Ez1(i))^2);
end
x=-3:0.1:3;
figure('name','Ex')
plot(x,Ex(:));
xlabel('x');
ylabel('Ex');
figure('name','Ey')
plot(x,Ey(:));
xlabel('x');
ylabel('Ey');
figure('name','Ez')
plot(x,Ez(:));
xlabel('x');
ylabel('Ez');


