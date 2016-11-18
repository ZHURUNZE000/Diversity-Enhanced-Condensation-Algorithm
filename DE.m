clc;
clear;
t0 = cputime;  
%差分?化算法程序  
%F0是??率 %Gm 最大迭代次数  
Gm = 3000;  
F0 = 0.5;  
Np = 100;  
CR = 0.9;  %交叉概率  
G= 1; %初始化代数  
D = 1; %所求??的?数  
Gmin = zeros(1,Gm); %各代的最??  
best_x = zeros(Gm,D); %各代的最?解  
value = zeros(1,Np);  
  
%?生初始?群  
%xmin = -10; xmax = 100;%??数的下界  
xmin = -5.12;  
xmax = 5.12;  


X0 = (xmax-xmin)*rand(Np,D) + xmin;  %?生Np个D?向量  
XG = X0;  
  
%%%%%%%%%%%%%----?里未做?价，不判断?止条件----%%%%%%%%%%%%%%%%%%%%%%%%  
  
XG_next_1= zeros(Np,D); %初始化  
XG_next_2 = zeros(Np,D);  
XG_next = zeros(Np,D);  
  
while G <= Gm   
G;   
%%%%%%%%%%%%%%%%%%%%%%%%----??操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = 1:Np  
        %?生j,k,p三个不同的数  
        a = 1;  
        b = Np;  
        dx = randperm(b-a+1) + a- 1;  
        j = dx(1);  
        k = dx(2);  
        p = dx(3);  
        %要保?与i不同  
        if j == i  
            j  = dx(4);  
            else if k == i  
                 k = dx(4);  
                else if p == i  
                    p = dx(4);  
                    end  
                end  
         end  
          
        %??算子  
        suanzi = exp(1-Gm/(Gm + 1-G));  
        F = F0*2.^suanzi;  
        %??的个体来自三个随机父代  
         
        son = XG(p,:) + F*(XG(j,:) - XG(k,:));         
        for j = 1: D  
            if son(1,j) >xmin  & son(1,j) < xmax %防止??超出?界  
                XG_next_1(i,j) = son(1,j);  
            else  
                XG_next_1(i,j) = (xmax - xmin)*rand(1) + xmin;  
            end  
        end  
    end  
   %%%%%%%%%%%%%%%%%%%%%%%---交叉操作----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
      
    for i = 1: Np  
        randx = randperm(D);% [1,2,3,...D]的随机序列     
        for j = 1: D  
              
            if rand > CR & randx(1) ~= j % CR = 0.9   
                XG_next_2(i,j) = XG(i,j);  
            else  
                XG_next_2(i,j) = XG_next_1(i,j);  
            end  
        end  
    end  
      
    %%%%%%%%%%%%%%%%%%----??操作---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = 1:Np  
        if f(XG_next_2(i,:)) < f(XG(i,:))  
              
            XG_next(i,:) = XG_next_2(i,:);  
        else  
            XG_next(i,:) = XG(i,:);  
        end  
    end  
      
    %找出最小?  
    for i = 1:Np  
        value(i) = f(XG_next(i,:));  
    end  
    [value_min,pos_min] = min(value);  
      
    %第G代中的目?函数的最小?  
    Gmin(G) = value_min;     
    %保存最?的个体  
    best_x(G,:) = XG_next(pos_min,:);     
      
    XG = XG_next;   %更新   ！！！
    trace(G,1) = G;  
    trace(G,2) = value_min;  
    G = G + 1;  
    
end  
  [value_min,pos_min] = min(Gmin);  
  best_value = value_min  
  best_vector =  best_x(pos_min,:)    
  fprintf('DE所耗的???：%f \n',cputime - t0);  
  %画出代数跟最?函数?之?的?系?    
  plot(trace(:,1),trace(:,2));  