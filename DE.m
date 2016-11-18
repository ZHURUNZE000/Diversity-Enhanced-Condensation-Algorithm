clc;
clear;
t0 = cputime;  
%����?���Z�@����  
%F0��??�� %Gm �ő�R�㎟��  
Gm = 3000;  
F0 = 0.5;  
Np = 100;  
CR = 0.9;  %�����T��  
G= 1; %���n���㐔  
D = 1; %����??�I?��  
Gmin = zeros(1,Gm); %�e��I��??  
best_x = zeros(Gm,D); %�e��I��?��  
value = zeros(1,Np);  
  
%?�����n?�Q  
%xmin = -10; xmax = 100;%??���I���E  
xmin = -5.12;  
xmax = 5.12;  


X0 = (xmax-xmin)*rand(Np,D) + xmin;  %?��Np��D?����  
XG = X0;  
  
%%%%%%%%%%%%%----?������?���C�s���f?�~����----%%%%%%%%%%%%%%%%%%%%%%%%  
  
XG_next_1= zeros(Np,D); %���n��  
XG_next_2 = zeros(Np,D);  
XG_next = zeros(Np,D);  
  
while G <= Gm   
G;   
%%%%%%%%%%%%%%%%%%%%%%%%----??����----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for i = 1:Np  
        %?��j,k,p�O���s���I��  
        a = 1;  
        b = Np;  
        dx = randperm(b-a+1) + a- 1;  
        j = dx(1);  
        k = dx(2);  
        p = dx(3);  
        %�v��?�^i�s��  
        if j == i  
            j  = dx(4);  
            else if k == i  
                 k = dx(4);  
                else if p == i  
                    p = dx(4);  
                    end  
                end  
         end  
          
        %??�Z�q  
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
    for i = 1:Np  
        value(i) = f(XG_next(i,:));  
    end  
    [value_min,pos_min] = min(value);  
      
    %��G�㒆�I��?�����I�ŏ�?  
    Gmin(G) = value_min;     
    %�ۑ���?�I����  
    best_x(G,:) = XG_next(pos_min,:);     
      
    XG = XG_next;   %�X�V   �I�I�I
    trace(G,1) = G;  
    trace(G,2) = value_min;  
    G = G + 1;  
    
end  
  [value_min,pos_min] = min(Gmin);  
  best_value = value_min  
  best_vector =  best_x(pos_min,:)    
  fprintf('DE���ՓI???�F%f \n',cputime - t0);  
  %��o�㐔���?����?�V?�I?�n?    
  plot(trace(:,1),trace(:,2));  