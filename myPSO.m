clc;clear;close all;
tic;                          %程序运行计时开始
E0=1.0e-20;                     %允许误差
MaxNum=100;                   %最大迭代次数
narvs=2;                      %目标函数的自变量个数，解空间维数
particlesize=100;              %粒子群规模
c1=2;                         %个体学习因子
c2=2;                         %社会学习因子
w=0.6;                        %惯性因子
vmax=1;                     %粒子的最大飞翔速度
rand('state',100);
x=-5 + 10*rand(particlesize,narvs); %粒子群的初始位置
v=0.5*rand(particlesize,narvs);     %粒子群的初始速度
%指定目标函数
fitness=inline('20+x.^2+y.^2-10*(cos(2*pi.*x)+cos(2*pi.*y))','x','y');      %
%fitness=inline('x.^2+y.^2','x','y');      

globalbest_faval_list = zeros(MaxNum,1);
%%%画出上述函数图像
[a,b]=meshgrid([-5:0.01:5]);
z = fitness(a,b);
figure();mesh(a,b,z)
xlabel('X1');ylabel('X2');title('函数值')
saveas(gcf, '函数值', 'png')
figure(2);contour(a,b,z);hold on 
xlabel('X1');ylabel('X2');title('等高线图')
saveas(gcf, '等高线图', 'png')
pic_num = 1;
%%%

for i=1:particlesize   
    f(i)=fitness(x(i,1),x(i,2));
end
personalbest_x=x;
personalbest_faval=f;
[globalbest_faval i]=min(personalbest_faval);
globalbest_x=personalbest_x(i,:);
k=1;
while k<=MaxNum
    %%%画粒子群位置图
    figure(3);
    scatter(x(:,1),x(:,2)) ; 
    axis([-5,5 -5,5])
    xlabel('X1');ylabel('X2');title(['Particles Location (iter=',num2str(k),')'])
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'particleLaction.gif','gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map,'particleLaction.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    pic_num = pic_num + 1;
    %%%更新位置和速度
    for i=1:particlesize
        f(i)=fitness(x(i,1),x(i,2));
        if f(i)<personalbest_faval(i)  %判断当前值是否是历史最佳
            personalbest_faval(i)=f(i);
            personalbest_x(i,:)=x(i,:);
        end
    end
    [globalbest_faval i]=min(personalbest_faval);
    globalbest_x=personalbest_x(i,:);
    for i=1:particlesize                    %更新每个个体的位置

        v(i,:)=w*v(i,:)+c1*rand*(personalbest_x(i,:)-x(i,:))+c2*rand*(globalbest_x-x(i,:));
        for j=1:narvs                       %判断是否超过最大速度
            if v(i,j)>vmax
                v(i,j)=vmax;
            elseif v(i,j)<-vmax
                v(i,j)=-vmax;
            end
        end
        x(i,:)=x(i,:)+v(i,:);
    end

    %%% globalbest_faval curve
    globalbest_faval_list(k,1) = globalbest_faval;
    
    if abs(globalbest_faval)<E0,break,end
    k=k+1;
end
globalbest_faval = num2str(globalbest_faval);
globalbest_x     = num2str(globalbest_x);
disp(strcat('最小值点','=',globalbest_x));
disp(strcat('最小值','=',globalbest_faval));
toc;


%globalbest_list plot
figure()
iter = 1:k-1;
plot(iter,globalbest_faval_list(1:k-1,:));
xlabel('iter');
ylabel('Globalbest Faval');
saveas(gcf, 'Globalbest Faval', 'png')




