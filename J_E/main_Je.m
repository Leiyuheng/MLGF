close all
clear
clc
%% ============================场源位置定义===========================
position_s = [0 0 -1.4];
position_f = [0 1 -0.3];
polar = [pi/6,7*pi/18,1]; % 这里的theta matlab函数中定义的是和xoy面的夹角

x_f = -3:0.1:3;
E = zeros(3,1,length(x_f));
for i = 1:length(x_f)
    position_f = [x_f(i) 1 -0.3];
    Ge = MLGF_Ge(position_s,position_f); % 场源相对位置（层数）确定 GF确定
    E(:,:,i) = calculate_E(Ge,polar);
end

Ex = squeeze(E(1, 1,:));
Ey = squeeze(E(2, 1,:));
Ez = squeeze(E(3, 1,:));
E_total = Ex.^2 + Ey.^2 +Ez.^2;

% 打开文件
fileID = fopen('J.efe', 'r');

% 读取文件头部信息，假设文件头部信息有固定行数
headerLines = 15; % 跳过行数
for i = 1:headerLines
    fgetl(fileID); % 逐行读取并忽略文件头部信息
end

% 读取数据
data = textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
% 关闭文件
fclose(fileID);

% 提取数据列
x_feko = data{1};
y_feko = data{2};
z_feko = data{3};

Ex_feko = data{4} + 1i*data{5};
Ey_feko = data{6} + 1i*data{7};
Ez_feko = data{8} + 1i*data{9};



% 绘制图形
figure;
hold on;

% 绘制仿真数据 (实线，绝对值)
plot(x_feko, abs(Ex_feko), 'r-', 'DisplayName', 'Ex (仿真)');
plot(x_feko, abs(Ey_feko), 'g-', 'DisplayName', 'Ey (仿真)');
plot(x_feko, abs(Ez_feko), 'b-', 'DisplayName', 'Ez (仿真)');

% 绘制计算数据 (点，绝对值)
plot(x_f, abs(Ex), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Ex (计算)');
plot(x_f, abs(Ey), 'g*', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Ey (计算)');
plot(x_f, abs(Ez), 'b+', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Ez (计算)');

% 设置图例和标签
legend show;
xlabel('Observation point x (m)');
ylabel('Magnetic field (A/m)');
title('Comparison of Simulated and Calculated Magnetic Fields');
grid on;
hold off;
