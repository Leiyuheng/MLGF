clc
clear
close all

position_f = [0 1 -0.3];
position_M1 = [1 0 -1.4];
position_M2 = [0 0 -1.4];
polar = [pi/6,7*pi/18,1];

x_f = linspace(-3,3,201);
H = zeros(3,1,length(x_f));
for i = 1:length(x_f)
    position_f = [x_f(i) 1 -0.3];
    Gh = MLGF_Gh(position_M1,position_f); % 场源相对位置（层数）确定 GF确定
    H(:,:,i) = calculate_H(Gh,polar);
    Gh = MLGF_Gh(position_M2,position_f); % 场源相对位置（层数）确定 GF确定
    H(:,:,i) = H(:,:,i) + calculate_H(Gh,polar);
end

Hx = squeeze(H(1, 1,:));
Hy = squeeze(H(2, 1,:));
Hz = squeeze(H(3, 1,:));
H_total = Hx.^2 + Hy.^2 +Hz.^2;


% 打开文件
fileID = fopen('2M.hfe', 'r');

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

Hx_feko = data{4} + 1i*data{5};
Hy_feko = data{6} + 1i*data{7};
Hz_feko = data{8} + 1i*data{9};



% 绘制图形
figure;
hold on;

plot(x_feko, abs(Hx_feko), 'r-', 'DisplayName', 'Hx (仿真)','LineWidth',1.5);
plot(x_feko, abs(Hy_feko), 'k-', 'DisplayName', 'Hy (仿真)','LineWidth',1.5);
plot(x_feko, abs(Hz_feko), 'b-', 'DisplayName', 'Hz (仿真)','LineWidth',1.5);

plot(x_f, abs(Hx), 'r*', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Hx (计算)');
plot(x_f, abs(Hy), 'k+', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'DisplayName', 'Hy (计算)');
plot(x_f, abs(Hz), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', 'Hz (计算)');

legend show;
xlabel('Observation point x (m)');
ylabel('Magnetic field (A/m)');
title('磁偶极子仿真与计算对比图');
grid on;
hold off;


error_Hx = abs(Hx-Hx_feko);
RMS_Hx = sqrt(mean(error_Hx.^2));
error_Hy = abs(Hy-Hy_feko);
RMS_Hy = sqrt(mean(error_Hy.^2));
error_Hz = abs(Hz-Hz_feko);
RMS_Hz = sqrt(mean(error_Hz.^2));
