%%This script is utilized to plot the 3D vectors of the field data.
% If the dataset is 3D, it will select a userdesigned page by the parameter
% zpage. The dataset need include the following parameters:
% Field: Ex1, Ey1, Ez1
% Axis: x_axis, y_axis, z_axis
% FileNum: used to tag data sources, therory or simulation(FDTD)
% 
% Coded by Sen Lu, Nanjing University of science and thehnology,
% lusen@njust.edu.cn
% Copyright (c) 2024, Sen Lu.
% Last update: 08/01/2025

clc
clear
close all

% load("BlochSkyrmionWithFocusedMethod_theory_cycl_center.mat");FileNum = 1;
load("BlochSkyrmionWithFocusedMethod_theory_cycl_right.mat");FileNum = 1;

% load("model_fs2.86_bw4_NAeff1.4_3D.mat");FileNum = 3;
% load("model_fs3.08_bw4_NAeff1.3_3D.mat");FileNum = 3;
% load("model_fs3.33_bw4_NAeff1.2_3D.mat");FileNum = 3;
% load("model_fs3.64_bw4_NAeff1.1_3D.mat");FileNum = 3;
% load("model_fs4.44_bw4_NAeff0.9_3D.mat");FileNum = 3;
% load("model_fs5_bw4_NAeff0.8_3D.mat");FileNum = 3;
% load("model_fs5.71_bw4_NAeff0.7_3D.mat");FileNum = 3;
% load("model_fs6.67_bw4_NAeff0.6_3D.mat");FileNum = 3;

%% 数据处理
% zpage = 3.1e-6; % um [2.3 2.7 2.9 3.1 3.5]
% pageNumber = find(abs(z_axis-zpage)==min(abs(z_axis-zpage)));

if FileNum == 3
    Ex1 = Ex1(:,:,pageNumber)';
    Ey1 = Ey1(:,:,pageNumber)';
    Ez1 = Ez1(:,:,pageNumber)';
end

% 坐标矩阵
xls = x_axis*1e6;
yls = y_axis*1e6;

% 单位化
I = real(Ex1).^2 + real(Ey1).^2 + real(Ez1).^2; % 单位化

reex = real(Ex1)./sqrt(I);
reey = real(Ey1)./sqrt(I);
reez = real(Ez1)./sqrt(I);

% 修正系数
cutcoef = 0.45;
axiscoef = 4.2/cutcoef;
%% 绘图

loop_Num = 7;                               % N是要画的箭头圈数
M = loop_Num-1;

for loop = 0:M
    loop_arrowNum = 3+4*loop;                  % 一圈有多少个箭头
    
    arrow_pos_rho = 1*loop*cutcoef/M;                   % 每一圈箭头对应的半径
    arrow_pos_theta = linspace(0,2*pi,loop_arrowNum);   % 每一个箭头对应的角度，角度等于2*pi/个数

    for loop_arrow_i = 1:loop_arrowNum                   % 绘制每一个箭头
        % 设置箭头起点位置
        [arrow_pos_x1, arrow_pos_y1] = pol2cart(arrow_pos_theta(loop_arrow_i), arrow_pos_rho);
        

        % 设置箭头终点位置
        % 求解箭头终点位置的索引，输出矩阵中某个值对应的x,y坐标（最小值法）
        dim1 = find(abs(xls-arrow_pos_x1)==min(abs(xls-arrow_pos_x1)), 1 );
        dim2 = find(abs(yls-arrow_pos_y1)==min(abs(yls-arrow_pos_y1)), 1 );

        arrow_pos_x2 = reex(dim2,dim1);
        arrow_pos_y2 = reey(dim2,dim1);
        arrow_pos_z2 = reez(dim2,dim1);

        % 设置箭头位置
        arrow_pos1 = [axiscoef*arrow_pos_x1,axiscoef*arrow_pos_y1,-0.5*arrow_pos_z2];
        arrow_pos2 = [arrow_pos_x2,arrow_pos_y2,arrow_pos_z2];

        % 绘制箭头
        arrow3D_up(arrow_pos1,arrow_pos2,arrow_pos_z2)     %第一个元素使箭头杆端的位置，第二个元素是箭头尖端相对杆端的向量，第三个元素是颜色值，这里取尖端的值来对应颜色
    end
    
end


% 颜色
axis equal tight off        % 坐标轴比例尺相同，且关掉了坐标轴显示
caxis([-1 1])                % 调整颜色区间
colormap("jet")               %给了一个color
view(0,37)

% 光照、材质
material dull
h1 = light;                  %给图像打一束光
lightangle(h1,50,0);
lighting gouraud;