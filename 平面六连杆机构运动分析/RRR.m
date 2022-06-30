%见教材(哈工大机械原理第七版)P37页
%输入
%Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy, -----B,D端点的位置，速度，加速度
%LBC,LCD----杆长
%Flag---标识符0，BCD顺时针排列，1--BCD逆时针排列
%输出
%Cx,Cy,Cdx,Cdy,Cddx,Cddy-----C端点的位置，速度，加速度
%theta_BC, theta_DC,omega_BC,omega_DC,alpha_BC,alpha_DC---两根杆角度、角速度和角加速度

function [Cx,Cy,Cdx,Cdy,Cddx,Cddy,theta_BC, theta_DC,omega_BC,omega_DC,alpha_BC,alpha_DC] = ...
RRR(Bx,By,Bdx,Bdy,Bddx,Bddy,Dx,Dy,Ddx,Ddy,Dddx,Dddy,LBC,LCD,Flag)

LBD = sqrt((Dx-Bx)*(Dx-Bx) + (Dy-By)*(Dy-By));

A = 2*LBC*(Dx-Bx);
B = 2*LBC*(Dy-By);
C = LBC*LBC + LBD*LBD - LCD*LCD;

%顺时针布置
if Flag == 0
    theta_BC = atan2(B,A) - atan2(-sqrt(A*A + B*B - C*C),C);
    %theta_BC = 2*atan( (B + sqrt(A*A + B*B -C*C))/(A + C));
end

%逆时针布置
if Flag == 1
    theta_BC = 2*atan( (B - sqrt(A*A + B*B -C*C))/(A + C));
end

Cx = Bx + LBC*cos(theta_BC);
Cy = By + LBC*sin(theta_BC);
theta_DC = atan2((Cy - Dy),(Cx - Dx));

%速度方程
CBC = LBC*cos(theta_BC);
SBC = LBC*sin(theta_BC);
CDC = LCD*cos(theta_DC);
SDC = LCD*sin(theta_DC);
G1 = CBC*SDC - CDC*SBC;

omega_BC = (CDC*(Ddx - Bdx) + SDC*(Ddy - Bdy))/G1;
omega_DC = (CBC*(Ddx - Bdx) + SBC*(Ddy - Bdy))/G1;

Cdx = Bdx - omega_BC*LBC*sin(theta_BC);
Cdy = Bdy + omega_BC*LBC*cos(theta_BC);

%加速度
G2 = Dddx - Bddx + omega_BC*omega_BC*CBC - omega_DC*omega_DC*CDC;
G3 = Dddy - Bddy + omega_BC*omega_BC*SBC - omega_DC*omega_DC*SDC;

alpha_BC = (G2*CDC + G3*SDC)/G1;
alpha_DC = (G2*CBC + G3*SBC)/G1;

Cddx = Bddx - omega_BC*omega_BC*LBC*cos(theta_BC) - alpha_BC*LBC*sin(theta_BC);
Cddy = Bddy - omega_BC*omega_BC*LBC*sin(theta_BC) + alpha_BC*LBC*cos(theta_BC);

end
