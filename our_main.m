
clc
clear all
%---------------导入决策矩阵------%
X=load('data2.txt');
w=[0.3 0.4 0.3];%属性权重
for dd=1:91
%     X(1,1)=(dd-1)/100;%-动态实验1----------------------%
   X(3,1)=(dd-1)/100;%-动态实验2----------------------%
%---------------计算加权决策矩阵------%
n=4;   %------方案的个数
k=5;   %------专家的个数
m=3;   %------属性指标的个数
for j=1:(n*k)
     for i=1:m
         Y(j,1+(i-1)*4)=1-(1-X(j,1+(i-1)*4))^w(i);
         Y(j,2+(i-1)*4)=1-(1-X(j,2+(i-1)*4))^w(i);
         Y(j,3+(i-1)*4)=(X(j,3+(i-1)*4))^w(i);
         Y(j,4+(i-1)*4)=(X(j,4+(i-1)*4))^w(i);
     end
end
WY=Y;
%-------------------------------------------------------%

%---------------计算正理想----------------%
for jj=1:k
        for ii=1:m
        Z(jj,1+4*(ii-1))=max(Y([jj jj+k jj+k*2 jj+k*3],1+4*(ii-1)));
        Z(jj,2+4*(ii-1))=max(Y([jj jj+k jj+k*2 jj+k*3],2+4*(ii-1)));
        Z(jj,3+4*(ii-1))=min(Y([jj jj+k jj+k*2 jj+k*3],3+4*(ii-1)));
        Z(jj,4+4*(ii-1))=min(Y([jj jj+k jj+k*2 jj+k*3],4+4*(ii-1)));
        end
end
ZZ=Z;

%---------------计算遗憾和满意矩阵----------------%
 for jjj=1:n
     A=Y(1+(jjj-1)*k:jjj*k,:);
         for j=1:k
             for i=1:m
                 %---------------计算遗憾矩阵----------------%
                 R(j+(jjj-1)*k,1+(i-1)*4)=Z(j,1+(i-1)*4)*A(j,3+(i-1)*4);
                 R(j+(jjj-1)*k,2+(i-1)*4)=Z(j,2+(i-1)*4)*A(j,4+(i-1)*4);
                 R(j+(jjj-1)*k,3+(i-1)*4)=Z(j,3+(i-1)*4)+A(j,1+(i-1)*4)-Z(j,3+(i-1)*4)*A(j,1+(i-1)*4);
                 R(j+(jjj-1)*k,4+(i-1)*4)=Z(j,4+(i-1)*4)+A(j,2+(i-1)*4)-Z(j,4+(i-1)*4)*A(j,2+(i-1)*4);
          
             end
         end
 end
 
%---------------计算正理想遗憾矩阵----------------%
 for jj=1:k
        for ii=1:m
        ZR(jj,1+4*(ii-1))=min(R([jj jj+k jj+k*2 jj+k*3],1+4*(ii-1)));
        ZR(jj,2+4*(ii-1))=min(R([jj jj+k jj+k*2 jj+k*3],2+4*(ii-1)));
        ZR(jj,3+4*(ii-1))=max(R([jj jj+k jj+k*2 jj+k*3],3+4*(ii-1)));
        ZR(jj,4+4*(ii-1))=max(R([jj jj+k jj+k*2 jj+k*3],4+4*(ii-1)));
        end
 end
%--------------------------------------------%
  %--------------计算正理想遗憾矩阵的模------------------------------%
 Z2=zeros(1,k);
 for kk=1:k
        for i=1:m
          Z2(1,kk)=Z2(1,kk)+(Z(kk,1+(i-1)*4))^2+(Z(kk,2+(i-1)*4))^2+(Z(kk,3+(i-1)*4))^2+(Z(kk,4+(i-1)*4))^2+(1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))^2+(1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))^2;
        end
 end

 Z22=sum(Z2);
 
%--------------------------------------------%
Y2=zeros(n,k);
YZ=zeros(n,k);
DYZ=zeros(n,k);
for j=1:n
    YY=Y(1+(j-1)*k:j*k,:);
    for kk=1:k
       for i=1:m
          Y2(j,kk)=Y2(j,kk)+(YY(kk,1+(i-1)*4))^2+(YY(kk,2+(i-1)*4))^2+(YY(kk,3+(i-1)*4))^2+(YY(kk,4+(i-1)*4))^2+(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4))^2+(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4))^2;
          YZ(j,kk)=YZ(j,kk)+(Z(kk,1+(i-1)*4))*(YY(kk,1+(i-1)*4))+(Z(kk,2+(i-1)*4))*(YY(kk,2+(i-1)*4))+(Z(kk,3+(i-1)*4))*(YY(kk,3+(i-1)*4))+(Z(kk,4+(i-1)*4))*(YY(kk,4+(i-1)*4))+(1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))*(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4))+(1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))*(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4));     
          DYZ(j,kk)=DYZ(j,kk)+(((Z(kk,1+(i-1)*4))-(YY(kk,1+(i-1)*4)))^2+((Z(kk,2+(i-1)*4))-(YY(kk,2+(i-1)*4)))^2+((Z(kk,3+(i-1)*4))-(YY(kk,3+(i-1)*4)))^2+((Z(kk,4+(i-1)*4))-(YY(kk,4+(i-1)*4)))^2+((1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))-(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4)))^2+((1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))-(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4)))^2)^0.5;
       end
    end
    Y22(j)=sum(Y2(j,:));
    YZZ(j)=sum(YZ(j,:));
    DDYZ(j)=sum(DYZ(j,:));
end
ZZZ=YZ;
%--------------------------------------------%   
 for j=1:n
    GR1(j)=YZZ(j)/Z22^0.5;
    GR2(j)=YZZ(j)/(YZZ(j)+abs(Y22(j)- Z22));
    GR3(j)=YZZ(j)/((Y22(j))^0.5*(Z22)^0.5+abs(Y22(j)- Z22)+abs(abs(YZZ(j)- Z22)));
 end
%--------------------------------------------%
 GN=[GR1;GR2;GR3;DDYZ];
 for j=1:n
    GCU1(j)=(GR1(j)-min(GR1))/(max(GR1)-min(GR1));
    GCU2(j)=(GR2(j)-min(GR2))/(max(GR2)-min(GR2));
    GCU3(j)=(GR3(j)-min(GR3))/(max(GR3)-min(GR3));
    GCU4(j)=(max(DDYZ)-DDYZ(j))/(max(DDYZ)-min(DDYZ));
    
 end
 
 
  Z=ZR;
  Y=[]; 
  Y=R;
  %--------------计算正理想遗憾矩阵的模------------------------------%
 Z2=zeros(1,k);
 for kk=1:k
        for i=1:m
          Z2(1,kk)=Z2(1,kk)+(Z(kk,1+(i-1)*4))^2+(Z(kk,2+(i-1)*4))^2+(Z(kk,3+(i-1)*4))^2+(Z(kk,4+(i-1)*4))^2+(1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))^2+(1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))^2;
        end
 end
 Z22=sum(Z2);
%--------------------------------------------%
Y2=zeros(n,k);
YZ=zeros(n,k);
DYZ=zeros(n,k);
  for j=1:n
    YY=Y(1+(j-1)*k:j*k,:);
    for kk=1:k
       for i=1:m
          Y2(j,kk)=Y2(j,kk)+(YY(kk,1+(i-1)*4))^2+(YY(kk,2+(i-1)*4))^2+(YY(kk,3+(i-1)*4))^2+(YY(kk,4+(i-1)*4))^2+(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4))^2+(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4))^2;
          YZ(j,kk)=YZ(j,kk)+(Z(kk,1+(i-1)*4))*(YY(kk,1+(i-1)*4))+(Z(kk,2+(i-1)*4))*(YY(kk,2+(i-1)*4))+(Z(kk,3+(i-1)*4))*(YY(kk,3+(i-1)*4))+(Z(kk,4+(i-1)*4))*(YY(kk,4+(i-1)*4))+(1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))*(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4))+(1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))*(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4));
          DYZ(j,kk)=DYZ(j,kk)+(((Z(kk,1+(i-1)*4))-(YY(kk,1+(i-1)*4)))^2+((Z(kk,2+(i-1)*4))-(YY(kk,2+(i-1)*4)))^2+((Z(kk,3+(i-1)*4))-(YY(kk,3+(i-1)*4)))^2+((Z(kk,4+(i-1)*4))-(YY(kk,4+(i-1)*4)))^2+((1-(Z(kk,1+(i-1)*4))-Z(kk,3+(i-1)*4))-(1-(YY(kk,1+(i-1)*4))-YY(kk,3+(i-1)*4)))^2+((1-(Z(kk,2+(i-1)*4))-Z(kk,4+(i-1)*4))-(1-(YY(kk,2+(i-1)*4))-YY(kk,4+(i-1)*4)))^2)^0.5;
       end
    end
    Y22(j)=sum(Y2(j,:));
    YZZ(j)=sum(YZ(j,:));
    DDYZ(j)=sum(DYZ(j,:));
end
    
%--------------------------------------------%   
 for j=1:n

    GR1(j)=YZZ(j)/Z22^0.5;
    GR2(j)=YZZ(j)/(YZZ(j)+abs(Y22(j)- Z22));
    GR3(j)=YZZ(j)/((Y22(j))^0.5*(Z22)^0.5+abs(Y22(j)- Z22)+abs(abs(YZZ(j)- Z22)));
    
 end
 GR=[GR1;GR2;GR3;DDYZ];
%--------------------------------------------%

 for j=1:n
    GRU1(j)=(GR1(j)-min(GR1))/(max(GR1)-min(GR1));
    GRU2(j)=(GR2(j)-min(GR2))/(max(GR2)-min(GR2));
    GRU3(j)=(GR3(j)-min(GR3))/(max(GR3)-min(GR3));
    GRU4(j)=(max(DDYZ)-DDYZ(j))/(max(DDYZ)-min(DDYZ));
    
 end
SS1=[GCU1;GRU1;GCU2;GRU2;GCU3;GRU3;GCU4;GRU4];
%--------------------------------------------%
        for j=1:n
          Q1(dd,j)=0.5*GCU1(j)+0.5*GRU1(j);
          Q2(dd,j)=0.5*GCU2(j)+0.5*GRU2(j);
          Q4(dd,j)=0.5*GCU3(j)+0.5*GRU3(j);
          Q3(dd,j)=0.5*GCU4(j)+0.5*GRU4(j);
        end
end
 Q1=Q1';
 Q2=Q2';
 Q3=Q3';
 Q4=Q4';
%--------------------------------------------%
a=0:90;
%--------------------------------------------%
%----------------------动态实验1图的绘制----------------------%
% figure(1)
% subplot(1,4,1)
% plot(a,Q1(1,:),'-g',a,Q1(2,:),'--r',a,Q1(3,:),'-b',a,Q1(4,:),'--k');
% xlabel('\mu^{l}_{11}=\theta/100'),ylabel('Classical projection measure Q-value');
% axis([0 90 -0.02 1.1])
% 
% 
% subplot(1,4,2)
% plot(a,Q2(1,:),'-g',a,Q2(2,:),'--r',a,Q2(3,:),'-b',a,Q2(4,:),'--k');
% xlabel('\mu^{l}_{11}=\theta/100'),ylabel('Q-value in reference[53]');
% axis([0 90 -0.02 1.1])
% subplot(1,4,3)
% plot(a,Q3(1,:),'-g',a,Q3(2,:),'--r',a,Q3(3,:),'-b',a,Q3(4,:),'--k');
% xlabel('\mu^{l}_{11}=\theta/100'),ylabel('Euclidean distance measure Q-value');
% axis([0 90 -0.02 1.1])
% 
% subplot(1,4,4)
% plot(a,Q4(1,:),'-g',a,Q4(2,:),'--r',a,Q4(3,:),'-b',a,Q4(4,:),'--k');
% xlabel('\mu^{l}_{11}=\theta/100'),ylabel('Q-value in this paper');
% legend('A_1','A_2','A_3','A_4');
% axis([0 90 -0.02 1.1])
%  
%----------------------动态实验2图的绘制----------------------%
figure(1)
subplot(1,4,1)
plot(a,Q1(1,:),'-g',a,Q1(2,:),'--r',a,Q1(3,:),'-b',a,Q1(4,:),'--k');
xlabel('\mu^{l}_{31}=\delta/100'),ylabel('Classical projection measure Q-value');
axis([0 90 -0.02 1.1])


subplot(1,4,2)
plot(a,Q2(1,:),'-g',a,Q2(2,:),'--r',a,Q2(3,:),'-b',a,Q2(4,:),'--k');
xlabel('\mu^{l}_{31}=\delta/100'),ylabel('Q-value in reference [53]');
axis([0 90 -0.02 1.1])

subplot(1,4,3)
plot(a,Q3(1,:),'-g',a,Q3(2,:),'--r',a,Q3(3,:),'-b',a,Q3(4,:),'--k');
xlabel('\mu^{l}_{31}=\delta/100'),ylabel('Euclidean distance measure Q-value');
axis([0 90 -0.02 1.1])

subplot(1,4,4)
plot(a,Q4(1,:),'-g',a,Q4(2,:),'--r',a,Q4(3,:),'-b',a,Q4(4,:),'--k');
xlabel('\mu^{l}_{31}=\delta/100'),ylabel('Q-value in this paper');
legend('A_1','A_2','A_3','A_4');
axis([0 90 -0.02 1.1])

