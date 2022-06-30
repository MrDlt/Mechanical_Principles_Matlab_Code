clc
clear
L1 = 26.5;
L2 = 105.6;
L21 = 65;
L3 = 67.5;
L4 = 87.5;
L5 = 34.4;
L6 = 25;
alpha = 60;


omega1 = 1; %rad/s
alpha1 = 0;
dr = pi/180; %ratio of deg. to rad.

Ax = 0; Ay = 0; Adx = 0; Ady = 0; Addx = 0; Addy = 0;
Dx = L4; Dy =0; Ddx = 0; Ddy = 0; Dddx = 0; Dddy = 0;
Gx = 153.5; Gy = 41.7; Gdx = 0; Gdy = 0; Gddx = 0; Gddy = 0;

phi = 0;
dphi =0;
ddphi = 0;

%rad. to deg.
rd = 180/pi;
deg = [0:1:360];
m = length(deg);

for n = 1:m
    theta1 = deg(n)*dr;
    [Bx,By,Bdx,Bdy,Bddx,Bddy] = RR(Ax,Ay,Adx,Ady,Addx,Addy,theta1,omega1,alpha1,L1);

    [Cx,Cy,Cdx,Cdy,Cddx,Cddy,theta_2(n), theta_3(n),omega_2, omega_3,alpha_2, alpha_3] = ...
    RRR(Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,L2,L3,0)
    theta_4(n) = theta_2(n)-60*dr;

    [Ex(n),Ey(n),Edx,Edy,Eddx,Eddy] = RR(Cx,Cy,Cdx,Cdy,Cddx,Cddy,theta_4(n),omega_2,alpha_2,L21);
    
    [Fx,Fy,Fdx,Fdy,Fddx,Fddy,theta_5(n), theta_6,omega_5, omega_6,alpha_5, alpha_6] = ...
    RRR(Ex(n),Ey(n),Edx,Edy,Eddx,Eddy,Gx,Gy,Gdx,Gdy,Gddx,Gddy,L5,L6,1);
    s(n) = sqrt((Ex(n)-Ex(1))*(Ex(n)-Ex(1))+(Ey(n)-Ey(1))*(Ey(n)-Ey(1)));
    v(n) = sqrt(Edx*Edx+Edy*Edy);
    a(n) = sqrt(Eddx*Eddx+Eddy*Eddy);
    
end

figure(1)
plot(deg,s, 'b');
legend('s');
title('刨刀位移线图');
xlabel('转角\theta_1/\circ')
ylabel('位移mm')
grid on; hold on;

figure(2)
plot(deg,v, 'r');
legend('v');
title('刨刀速度线图');
xlabel('原动件转角\theta_1/\circ')
ylabel('速度mm/s')
grid on; hold on;

figure(3)
plot(deg,a, 'b');
legend('a');
title('加速度线图');
xlabel('原动件转角\theta_1/\circ')
ylabel('加速度mm/s^2')
grid on; hold on;

figure(4)
plot(Ex,Ey, 'r');
legend('Ey');
title('E点轨迹图');
xlabel('Ex mm')
ylabel('Ey mm')
grid on; hold on;

figure(5)
mm = moviein(30);
j = 0;
for n = 1:m
    j = j + 1;
    clf;

    x(1) = 0;
    y(1) = 0;
    x(2) = L1*cos(deg(n));
    y(2) = L1*sin(deg(n));
    x(3) = L4 + L3*cos(deg(n));
    y(3) = L3*sin(deg(n));
    x(4) = L4;
    y(4) = 0;
    x(5) = x(3) + L21*cos(theta_4(n));
    y(5) = y(3) + L21*sin(theta_4(n));
    x(6) = x(5) + L5*cos(theta_5(n));
    y(6) = y(5) + L5*sin(theta_5(n));
    x(7) = Gx;
    y(7) = Gy;


    plot(x,y);
    grid on; 
    plot(x(1),y(1),'o');
    hold on;
    plot(x(2),y(2),'o');
    hold on;
    plot([x(1),x(2),x(3),x(5),x(6),x(7)],[y(1),y(2),y(3),y(5),y(6),y(7)]);
    hold on;
    plot(x(3),y(3),'o');
    hold on;
    plot(x(4),y(4),'o');
    hold on;
    plot([x(1),x(4),x(3)],[y(1),y(4),y(3)]);
    hold on;
    plot(x(5),y(5),'o');
    hold on;     
    plot(x(6),y(6),'o');
    hold on;
    plot(x(7),y(7),'o');
    hold on;

    axis([-150 350 -150 200]);
    
    mm(j) = getframe;
end
movie(mm,3);