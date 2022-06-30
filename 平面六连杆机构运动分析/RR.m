%见教材(哈工大机械原理第七版)P36
%输入
%x,y,dx,dy,ddx,ddy -----A端点的位置，速度，加速度
%theta,omega,alpha----杆件的角度，角速度，角加速度theta,omega,alpha
%len----杆长
%输出
%Bx,By,Bdx,Bdy,Bddx,Bddy-----B端点的位置，速度，加速度

function [Bx,By,Bdx,Bdy,Bddx,Bddy] = RR(Ax,Ay,Adx,Ady,Addx,Addy,theta,omega,alpha,Len)
%位置
Bx = Ax + Len*cos(theta);
By = Ay + Len*sin(theta);

%速度
Bdx = Adx - omega*Len*sin(theta);
Bdy = Ady + omega*Len*cos(theta);

%加速度
Bddx = Addx - omega*omega*Len*cos(theta) - alpha*Len*sin(theta);
Bddy = Addy - omega*omega*Len*sin(theta) + alpha*Len*cos(theta);

end
