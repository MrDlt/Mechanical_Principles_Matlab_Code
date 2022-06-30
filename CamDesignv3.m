%偏置直动滚子推杆盘形凸轮轮廓曲线设计
%初始基圆半径20mm,滚子半径14mm,推程角165，远休角50，回程运动角100，近休角45.
%推杆推程等加速等减速运动，行程30mm,回程以余弦加速运动返回原处
%逆时针回转，凸轮回转中心偏于推杆左侧
%推程压力角30度，回程压力角75度（设计时，要调整基圆半径，在推程，回程时满足许应的压力角）。
%偏距e与基圆半径的比例为1:2

%start
clc; clear;

rd = 180/pi; %rad -> deg弧度与角度转换系数
dr = pi/180;

%已知参数
r0 = 20; %基圆半径，初始值
rr = 14; %滚子半径
h = 30; %行程
e = 1/2*r0; %偏距
deltar0 = 1; %基圆半径增加值
alpha1allow = 30*dr;
alpha2allow = 75*dr;
rhominallow = 0.3*rr;

delta1 = 165*dr; %转换成弧度
delta2 = 50*dr;
delta3 = 100*dr;
delta4 = 45*dr;

%delta12 表示两部分累加的角度,其他类似
delta12 = delta1 + delta2;
delta13 = delta1 + delta2 + delta3;
delta14 = delta1 + delta2 + delta3 + delta4;

%假定凸轮按照间隔角度转动
deltaDeg = 10; %间隔10度
omega1 = 1; %凸轮转速
deg = [0:deltaDeg:360]; %凸轮转动度数
N = length(deg); %转动1周，计算的点数

% %程序中变量的说明
% rhomin = 1000;deltarhomin = 0; %理论轮廓最小曲率半径及对用的凸轮转角
% alpha1max = 0; deltaalpha1max = 0; % 推程最最大压力角及对用的凸轮转角弧度
% alpha2max = 0; deltaa2pha1max = 0; % 回程最最大压力角及对用的凸轮转角弧度
% %x y xr xr 分别为理论和实际轮廓数据点

%从动件运动规律
while 1
    e = 1/2*r0;
    s0 = sqrt(r0*r0 - e*e);

    rhomin = 1000;deltarhomin = 0;
    alpha1max = 0; deltaalpha1max = 0;
    alpha2max = 0; deltaa2pha1max = 0;

    for n = 1:N
        %推程阶段
        rdeg = deg(n)*dr; %转动角度转化为弧度

        
        if rdeg <= delta1/2
            s(n) = 2*h*(rdeg/delta1)*(rdeg/delta1);
            %在这里，v不表示推杆的速度，而是表示类速度ds/ddelta，要注意。
            v(n) = 4*h*rdeg/delta1/delta1; %旧版本这里有问题
            ds = v(n);
            %计算推程压力角
            alpha1 = atan((abs(ds) - e)/(s(n) + s0));
            %选出推程最大的压力角
            if alpha1 > alpha1max
                alpha1max = alpha1;
                deltaalpha1max = rdeg;
            end  

        elseif (rdeg > delta1/2) && (rdeg <= delta1)
            s(n) = h-2*h/delta1/delta1*(delta1-rdeg)*(delta1-rdeg);
            %在这里，v不表示推杆的速度，而是表示类速度ds/ddelta，要注意。
            v(n) = 4*h*(delta1-rdeg)/delta1/delta1; %旧版本这里有问题
            ds = v(n);
            %计算推程压力角
            alpha1 = atan((abs(ds) - e)/(s(n) + s0));
            %选出推程最大的压力角
            if alpha1 > alpha1max
                alpha1max = alpha1;
                deltaalpha1max = rdeg;
            end  
                  
        %远休角
        elseif rdeg > delta1 && rdeg <= delta12
            s(n) = h; v(n) = 0;
            ds = v(n);
        %回程
        elseif rdeg > delta12 && rdeg <= delta13
            degback = rdeg - delta12;
            s(n) = 0.5*h*(1 + cos(pi*degback/delta3));
            v(n) = -0.5*pi*h*sin(pi*degback/delta3)/(delta3); %旧版本这里有问题
            ds = v(n);
            %计算回程压力角
            alpha2 = atan((abs(ds) + e)/(s(n) + s0));
            if alpha2 > alpha2max
                alpha2max = alpha2;
                deltaalpha2max = rdeg;
            end
        %近休角
        elseif rdeg > delta13 && rdeg <= delta14
            s(n) = 0; v(n) = 0;
            ds = v(n);
        end
    
        rho = 0; %公式
        if rho < rhomin
            rhomin = rho;
            deltarhomin = rdeg;
        end
        %------计算凸轮轮廓曲线------------------------------------------
        %计算理论轮廓曲线
        x(n) = (s0 + s(n))*sin(rdeg) + e*cos(rdeg);
        y(n) = (s0 + s(n))*cos(rdeg) - e*sin(rdeg);
        %对delta的导数
        dx(n) = (ds - e)*sin(rdeg) + (s0 + s(n))*cos(rdeg);
        dy(n) = (ds - e)*cos(rdeg) - (s0 + s(n))*sin(rdeg);
        %计算实际轮廓曲线
        stheta = dx(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        ctheta = -dy(n)/(sqrt(dx(n)*dx(n) + dy(n)*dy(n)));
        %内包络轮廓,用“-”号 %旧版本这里有问题（注释有误）
        xr(n) = x(n) - rr*ctheta;
        yr(n) = y(n) - rr*stheta;

    end %for

    %如果不满足设计参数，可以调整基圆半径的大小
    if alpha1max > alpha1allow || alpha2max > alpha2allow
        r0 = r0 + deltar0;
        continue;
    else
        break;
    end
end %while

%打印相关参数参数，并画出凸轮轮廓图
fprintf('基圆半径\n');
fprintf('%6.4f\n', r0);
fprintf('推程最大压力角，相应凸轮转角\n');
fprintf('%6.4f %6.4f\n',alpha1max*rd,deltaalpha1max*rd);
fprintf('回程最大压力角，相应凸轮转角\n');
fprintf('%6.4f %6.4f\n',alpha2max*rd,deltaalpha2max*rd);

%--------输出理论轮廓数据------------------------------------------------
fprintf('Results: nominal profile points \n');
fprintf('n	    x	     y \n');
for i = 1:N
    fprintf('%d\t %6.4f\t %6.4f \n', i,x(i),y(i));
end

%--------输出实际轮廓数据------------------------------------------------
fprintf('Results: actual profile points \n');
fprintf('n	    x	     y \n');
for i = 1:N
    fprintf('%d\t %6.4f\t %6.4f \n', i,xr(i),yr(i));
end

%-------------画出理论轮廓线-----------------------------------------------
figure(1)
hold on; grid on; axis equal;
title('偏置直动滚子推杆盘形凸轮设计')
xlabel('x/mm');
ylabel('y/mm');
plot(x,y, 'r-');
ct = linspace(0, 2*pi);
plot(r0*cos(ct), r0*sin(ct), 'g-'); %基圆
plot(e*cos(ct), e*sin(ct), 'c-'); %偏置圆
%------------画出实际轮廓线------------------------------------------------

plot(xr,yr, 'b-');
%------------画出运动规律图------------------------------------------------
figure(2)
plot(deg,s,'r-')
%--------------------------------------------------------------------------