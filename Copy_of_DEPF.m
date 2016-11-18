clc;
clear all;
close all;
x = 0; %���n
R = 1;
Q = 1;
tf = 100; %���H
N = 300;  %���q����
P = 2;
xhatPart = x;
for i = 1 : N    
    xpart(i) = x + sqrt(P) * randn;%���n��?����0�ϕ��� sqrt(P)�I���z���z
end


Gm = 300;  
F0 = 0.5;  
Np = 100;  
CR = 0.9;  %�����T��  
G= 1; %���n���㐔  
D = 1; %�����I��  
%Gmin = zeros(1,Gm);   
%best_x = zeros(Gm,D); %�e��I��?��  
value = zeros(1,Np);  



xArr = [x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];


for k = 1 : tf    
    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
    %
    %x=xArr1(k);
    y = x^2 / 20 + sqrt(R) * randn;  %k?��???
    %y=yArr1(k);
    
 for i = 1 : N
     
     xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
         + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;%N�����q
     ypart = xpartminus(i)^2 / 20;%
     vhat = y - ypart;%
     q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); 
     %���q�I���R�������x
 end
 qsum = sum(q);
for i = 1 : N
    q(i) = q(i) / qsum;%???�ꉻ
end  

  for i = 1 : N %�d�V��
      u = rand;
      qtempsum = 0;
      for j = 1 : N
          qtempsum = qtempsum + q(j);
          if qtempsum >= u
              xpart(i) = xpartminus(j);
              break;
          end
      end
  end
xhatPart = mean(xpart);
%�ō@�I��?��??��?N�����q�I����?�C?��??�d�V��?�@�e�����q�I??����  
  
XG(:,1) =  xpart;
xmin = min(xpart);  
xmax = max(xpart);    
XG_next_1= zeros(N,D); %���n��  
XG_next_2 = zeros(N,D);  
XG_next = zeros(N,D);  
  
while G <= Gm   
%G;   
%%%%%%%%%%%%%%%%%%%%%%%%----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = 1:Np  
        %j,k,p�O���s���I��  
        a = 1;  
        b = Np;  
        dx = randperm(b-a+1) + a- 1;  
        j = dx(1);  
        k = dx(2);  
        p = dx(3);  
        %�^i�s��  
        if j == i  
            j  = dx(4);  
            else if k == i  
                 k = dx(4);  
                else if p == i  
                    p = dx(4);  
                    end  
                end  
         end  
          
        %�Z�q  
        suanzi = exp(1-Gm/(Gm + 1-G));  
        F = F0*2.^suanzi;  
        %??�I���̗����O����������  
         
        son = XG(p,:) + F*(XG(j,:) - XG(k,:));         
        for j = 1: D  
            if son(1,j) >xmin  & son(1,j) < xmax %�h�~??���o?�E  
                XG_next_1(i,j) = son(1,j);  
            else  
                XG_next_1(i,j) = (xmax - xmin)*rand(1) + xmin;  
            end  
        end  
    end  
   %%%%%%%%%%%%%%%%%%%%%%%---��������----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
      
    for i = 1: Np  
        randx = randperm(D);% [1,2,3,...D]�I��������     
        for j = 1: D  
              
            if rand > CR & randx(1) ~= j % CR = 0.9   
                XG_next_2(i,j) = XG(i,j);  
            else  
                XG_next_2(i,j) = XG_next_1(i,j);  
            end  
        end  
    end  
      
    %%%%%%%%%%%%%%%%%%----??����---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = 1:Np  
        if f(XG_next_2(i,:)) < f(XG(i,:))  
              
            XG_next(i,:) = XG_next_2(i,:);  
        else  
            XG_next(i,:) = XG(i,:);  
        end  
    end  
      
    %�Q�o�ŏ�?  
%     for i = 1:Np  
%         value(i) = f(XG_next(i,:));  
%     end  
%     [value_min,pos_min] = min(value);  
      
    %��G�㒆�I��?�����I�ŏ�?  
   % Gmin(G) = value_min;     
    %�ۑ���?�I����  
   % best_x(G,:) = XG_next(pos_min,:);     
      
    XG = XG_next;   %�X�V   �I�I�I
    %trace(G,1) = G;  
    %trace(G,2) = value_min;  
    G = G + 1;  
    
end  
  clear xpart;
  xpart=XG;
  
  
  
  

xArr = [xArr x];   
yArr = [yArr y];  
% xhatArr = [xhatArr xhat]; 
PArr = [PArr P]; 
xhatPartArr = [xhatPartArr xhatPart];
end
t = 0 : tf;
figure;
plot(t, xArr, 'b-.', t, xhatPartArr, 'k-');
legend('Real Value','Estimated Value');
set(gca,'FontSize',10); 
xlabel('time step'); 
ylabel('state');
title('DE and Particle filter')
xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
figure;
plot(t,abs(xArr-xhatPartArr),'b');
title('The error of DEPF')