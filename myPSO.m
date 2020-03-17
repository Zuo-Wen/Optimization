clc;clear;close all;
tic;                          %�������м�ʱ��ʼ
E0=1.0e-20;                     %�������
MaxNum=100;                   %����������
narvs=2;                      %Ŀ�꺯�����Ա�����������ռ�ά��
particlesize=100;              %����Ⱥ��ģ
c1=2;                         %����ѧϰ����
c2=2;                         %���ѧϰ����
w=0.6;                        %��������
vmax=1;                     %���ӵ��������ٶ�
rand('state',100);
x=-5 + 10*rand(particlesize,narvs); %����Ⱥ�ĳ�ʼλ��
v=0.5*rand(particlesize,narvs);     %����Ⱥ�ĳ�ʼ�ٶ�
%ָ��Ŀ�꺯��
fitness=inline('20+x.^2+y.^2-10*(cos(2*pi.*x)+cos(2*pi.*y))','x','y');      %
%fitness=inline('x.^2+y.^2','x','y');      

globalbest_faval_list = zeros(MaxNum,1);
%%%������������ͼ��
[a,b]=meshgrid([-5:0.01:5]);
z = fitness(a,b);
figure();mesh(a,b,z)
xlabel('X1');ylabel('X2');title('����ֵ')
saveas(gcf, '����ֵ', 'png')
figure(2);contour(a,b,z);hold on 
xlabel('X1');ylabel('X2');title('�ȸ���ͼ')
saveas(gcf, '�ȸ���ͼ', 'png')
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
    %%%������Ⱥλ��ͼ
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
    %%%����λ�ú��ٶ�
    for i=1:particlesize
        f(i)=fitness(x(i,1),x(i,2));
        if f(i)<personalbest_faval(i)  %�жϵ�ǰֵ�Ƿ�����ʷ���
            personalbest_faval(i)=f(i);
            personalbest_x(i,:)=x(i,:);
        end
    end
    [globalbest_faval i]=min(personalbest_faval);
    globalbest_x=personalbest_x(i,:);
    for i=1:particlesize                    %����ÿ�������λ��

        v(i,:)=w*v(i,:)+c1*rand*(personalbest_x(i,:)-x(i,:))+c2*rand*(globalbest_x-x(i,:));
        for j=1:narvs                       %�ж��Ƿ񳬹�����ٶ�
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
disp(strcat('��Сֵ��','=',globalbest_x));
disp(strcat('��Сֵ','=',globalbest_faval));
toc;


%globalbest_list plot
figure()
iter = 1:k-1;
plot(iter,globalbest_faval_list(1:k-1,:));
xlabel('iter');
ylabel('Globalbest Faval');
saveas(gcf, 'Globalbest Faval', 'png')




