%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              2-mass-2-spring system simulation program                  %
%                     Written by Guo Yiming                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  |                                      %
%                               /=====                                    %
%                               \=====\   <===k2                          %
%                                =====/                                   %
%                                  |                                      %
%                            ______________                               %
%                            |             |                              % 
%                  _____     |    m1       |                              %
%                    |       |             |                              %
%              u2    |       |_____________|                              %
%                    V              |                                     %
%                               /=====                                    %
%                               \=====\   <==k1                           %
%                                =====/                                   %
%                                  |                                      %
%                            ______________                               %
%                            |             |                              % 
%                  _____     |    m2       |                              %
%                   |        |             |                              %
%              u1   |        |_____________|                              %
%                   V                                                     %
%%
% *Introduction*
% This Program is written to solve the 2-mass-2-spring system with MATLAB
% simulation.The initial question is from the text book(Boyce,Elementary 
% differential equations and boundary value problems).Page 181 Question 29.
% The figure is attached here 
%        <<https://i.loli.net/2021/06/28/BafGyZwjikKlMoY.png>>
%If you need more information please
%% Prerequest to Run this Program(Symbolic Toolbox and Optimization Toolbox)
%Since the process uses the solution of multivariable equations,
%multivariable equation and other solver to assist.This program uses the
%*Symbolic Toolbox* and *Optimization Toolbox*. This two Tool box may not 
% Preinstall in you MATLAB. So if you want to run it, don't forget to
% install them before. 
%
%% Variable List
% m1     = mass of the top block in the graph
% m2     = mass of the buttom bock in the graph
% k1     = the spring coefficient of the upper spring
% k2     = the spring coefficient of the downer spring
% u10    = the initial displacement of m_1=
% u1d0    = The initial velocity of m_1=
% u20     = The initial displacement of m_2=
% u2d0    = The initial velocity of m_2=
%% Input the initial condition and the constant variable
% Please enter the values of the variables requested by the computer 
% in turn

clc
syms r c11 c12 c13 c14 t
m1=input('m_1=');
m2=input('m_2=');
k1=input('k_1=');
k2=input('k_2=');
u10=input('The initial displacement of m_1=');
u1d0=input('The initial velocity of m_1=');
u20=input('The initial displacement of m_2=');
u2d0=input('The initial velocity of m_2=');

%% Calculation Process-gerneate more initial velocity
% The second and third derivatives of u1,u2 with respect to time are 
% obtained by recursive equations and used to determine the specific 
% solution later.

u1dd0=-(k1+k2)/m1*u10+k2*u20; %The second derivative of u1
u2dd0=-k2/m2*u20+k2/m2*u10; %The second derivative of u2
u1ddd0=-(k1+k2)/m1*u1d0+k2*u2d0; %The third derivative of u1
u2ddd0=-k2/m2*u2d0+k2/m2*u1d0; %The third derivative of u2

%% Solve the coeficient of FSS
% Use the solve() function to find the roots of the 4th degree equation.
eqn=(m1*m2)/k2*r^4+(m2/k2*(k1+k2)+m1)*r^2+k1==0;
Y=solve(eqn,r);
%% Find the corresponding fundamental solution based on the root form
% If the root ri is a real number, then its corresponding FS is e^{ri}
% If the root contains an imaginary part, and the imaginary part is 
% positive, then the corresponding FS is e^{Re(ri)}*sin(abs(Im(ri))*t)
% If the root contains an imaginary part, and the imaginary part is 
% negative, then the corresponding FS is e^{Re(ri)}*cos(abs(Im(ri))*t)
if isreal(Y(1))
        FSS1=@(t) exp(double(Y(1))*t);
    else
        if sign(imag(Y(1))) == 1
            FSS1=exp(double(real(Y(1)))*t)*sin(abs(double(imag(Y(1))))*t);
        else FSS1=exp(double(real(Y(1)))*t)*cos(abs(double(imag(Y(1))))*t);
        end
end

if isreal(Y(4))
        FSS4=@(t) exp(double(Y(4))*t);
    else
        if sign(imag(Y(4))) == 1
            FSS4=exp(double(real(Y(4)))*t)*sin(abs(double(imag(Y(4))))*t);
        else FSS4=exp(double(real(Y(4)))*t)*cos(abs(double(imag(Y(4))))*t);
        end
end

if isreal(Y(2))
        FSS2=@(t) exp(double(Y(2))*t);
    else
        if sign(imag(Y(2))) == 1
            FSS2=exp(double(real(Y(2)))*t)*sin(abs(double(imag(Y(2))))*t);
        else FSS2=exp(double(real(Y(2)))*t)*cos(abs(double(imag(Y(2))))*t);
        end
end


if isreal(Y(3))
        FSS3=@(t) exp(double(Y(3))*t);
    else
        if sign(imag(Y(3))) == 1
            FSS3=exp(double(real(Y(3)))*t)*sin(abs(double(imag(Y(3))))*t);
        else FSS3=exp(double(real(Y(3)))*t)*cos(abs(double(imag(Y(3))))*t);
        end
end


%% Combine and generate a general homogenous solution
% Combine and generate a general homogenous solution and find the 
% derivative of the function

u1=c11*FSS1+c12*FSS2+c13*FSS3+c14*FSS4;
u1d=diff(u1,t);
u1dd=diff(diff(u1,t),t);
u1ddd=diff(diff(diff(u1,t)),t);


%% Get the u1 function
% Solve the quadratic equation, obtain the four coefficients 
% and determine the function of u1
eqns=[subs(u1,t,0)==u10,subs(u1d,t,0)==u1d0,subs(u1dd,t,0)==u1dd0,subs(u1ddd,t,0)==u1ddd0];
sss=solve(eqns ,[c11 c12 c13 c14]);

u1=double(sss.c11)*FSS1+double(sss.c12)*FSS2+double(sss.c13)*FSS3+double(sss.c14)*FSS4;
u1f=matlabFunction(u1);
%% Get the u2 function
%By using recursive equation to solve the u2 function.
u2=m1/k2*diff(diff(u1,t),t)+(k1+k2)/k2*u1;
u2f=matlabFunction(u2);
%% Plot the graph
%Sample the two function lines and draw the function image
t=[0:0.01:30];
u1fs=u1f(t);
u2fs=u2f(t);
plot(t,u1fs+5)
hold on
plot(t,u2fs)
hold on

%% Draw balance point auxiliary line
%
o1=5+0*t;
o2=0*t;
plot(t,o1,'--');
hold on
plot(t,o2,'--');
legend('block1','block2','equilibrium point of block1','equilibrium point of block2')
title(['$m_1=$',num2str(m1),'$\quad m_2=$',num2str(m2),'$\quad k_1=$',num2str(k1),'$\quad k_2=$',num2str(k2),'$\quad u_1(0)=$',num2str(u10),'$\quad u_2(0)=$',num2str(u20),'$\quad u_1^\prime(0)=$',num2str(u1d0),'$\quad u_2^\prime(0)=$',num2str(u2d0)],'interpreter','latex');
xlabel('t(s)')
ylabel('u(m)')



