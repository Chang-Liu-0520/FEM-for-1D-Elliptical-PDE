%This program is used to solve one-dimension elliptic equation by Finite Element Method
%Chang Liu
%2016/6/1
clear all;
close all;
clc;
N = 100; %This is the number of all elements
xleft = 0;
xright = 1;
x = linspace(xleft,xright,N+1); %initialize the x coordinate
h = diff(x);
a = zeros(2,2,N); %initialize the local stiffness matrix
f = zeros(2,1,N); %initialize the local force
A = zeros(N+1,N+1); %initialize the global stiffness matrix
F = zeros(N+1,1) %initialize the global force

for k =1:N;
	a(1,1,k) = (1+(x(k)^2+x(k)*x(k+1)+x(k+1)^2)/3)/h(k)+h(k)/3;
	a(1,2,k) = -(1+(x(k)^2+x(k)*x(k+1)+x(k+1)^2)/3)/h(k)+h(k)/6;
	a(2,1,k) = a(1,2,k);
	a(2,2,k) = a(1,1,k); %build the local stiffness matrix
	f(1,1,k) = 1/h(k)*x(k+1)*(-2*(x(k+1)-x(k))-5/3*(x(k+1)^3-x(k)^3))+1/h(k)*(x(k+1)^2-x(k)^2)+1/h(k)*5/4*(x(k+1)^4-x(k)^4);
	f(2,1,k) = -x(k)/h(k)*(-2*(x(k+1)-x(k))-5/3*(x(k+1)^3-x(k)^3))-1/h(k)*(x(k+1)^2-x(k)^2)-1/h(k)*5/4*(x(k+1)^4-x(k)^4);
end

a(2,2,N) = a(2,2,N)+1;
f(2,1,N)=f(2,1,N)+5;

for i = 1:N %obtain the global stiffness matrix and global force using local stiffness and force matrix
	A(i,i)=A(i,i)+a(1,1,i);
	A(i,i+1)=A(i,i+1)+a(1,2,i);
	A(i+1,i)=A(i+1,i)+a(2,1,i);
	A(i+1,i+1)=A(i+1,i+1)+a(2,2,i);
	F(i,1)=F(i,1)+f(1,1,i);
	F(i+1,1)=F(i+1,1)+f(2,1,i);
end

A_con = A(2:N+1,2:N+1);
F_con = F(2:N+1);

u_con=inv(A_con)*F_con; % get the coefficieint of before the shpe function
u_con_total = zeros(N+1,1);
u_con_total(2:N+1,1)=u_con;
x_con = x(2:N+1);

n=10001
x_plot = linspace(0,1,n);
u_plot = zeros(1,n); %the u value calculated by FEM
y_plot = x_plot.^2; %the real value according to analysis solve

for i = 2:n
	left = max(find(x_plot(i)>x));
	right = left+1;
	hk = x(right)-x(left);
	u_plot(i)=u_con_total(left)*(x(right)-x_plot(i))/hk+u_con_total(right)*(x_plot(i)-x(left))/hk;
end
figure(1);
plot(x_plot,y_plot,'linewidth',2);
ylabel('u','FontWeight','bold');
xlabel('x','FontWeight','bold');
hold on;
plot(x_plot,u_plot,'linewidth',2)
hold on;
plot(x_con,u_con,'o','MarkerSize',8);
title('The Comparision between FEM solution and analytical solution');
legend('The analytical solution: u(x)=x^2', 'The numerical solution by FEM u_h(x)','The value at each lattice point');

error = y_plot-u_plot;
figure(2);
plot(x_plot,error);
ylabel('Absolute error','FontWeight','bold');
xlabel('x','FontWeight','bold');
%title('The Absolute error between FEM solution and analytical solution','position',[400,1000]);
error2 = sqrt(dot(error,error));

error_each=[0.1702,0.0017,1.7003e-5];
number = [1/10,1/100,1/1000];
figure(3);
plot(-log10(error_each),-log10(number));
hold on;
plot(-log10(error_each),-log10(number),'o');
xlabel('-log_{10}(h)','FontWeight','bold');
ylabel('-log_{10}(||u-u_h||)','FontWeight','bold');
title('The correlationship between ||u-u_h|| and space step h');