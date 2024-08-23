clear;
clc;
%%参数初始化
h_max=1;
h_min=0;
ub=100;
lb=-100;
E0=100;
c=0.2;
beta=1.5;
T_max=5000;
n=50;
d=30;
delta=1.42;
X=rand(n,d,T_max);
fit_best=zeros(1,T_max);
fit_best_global=zeros(1,T_max);
X_best_global=rand(T_max,d);
X_best=rand(T_max,d);
index=rand(1,T_max);
index_rank=zeros(1,5);
X_best_rank=zeros(5,d);
fit_rank=zeros(1,n);
funclabel=10
%%程序主体
%X和H的初始化

for i=1:n
 
 R_0=unifrnd(-1,1,1,d);
 Tar(i,:,1)=(1/(sqrt(dot(R_0,R_0))))*R_0;
 X(i,:,1)=lb*ones(1,d)+rand(1,d)*( ub-lb);
 H(i,:,1)=(h_min+rand*(h_max-h_min))*Tar(i,:,1);
 
 fit(1,i)=cec17_func((X(i,:,1))',funclabel);
end 
%寻找局部和全局适应度以及最优解
[fit_best(1,1),index(1,1)]=min(fit(1,:));
fit_best_global(1,1)=fit_best(1,1);
X_best(1,:)=X(index(1,1),:,1);
X_best_global(1,:)=X_best(1,:);
fit_rank(1,:)=sort(fit(1,:));
r5=0.49+rand*(1-0.49);
for t=1:T_max-1
  gravity(1,:)=unifrnd(-1,1,1,d);%生成重力向量
  Tar_1(1,:)=unifrnd(-1,1,1,d);
  H_target(1,:)=(1/(sqrt(dot(Tar_1,Tar_1))))*Tar_1;%生成目标磁场   
  lamda=-E0*(1-(power(t,0.48)/T_max))*log(1-(power(t,0.48)/T_max));
  lamda_find=lamda;
    for i=1:25
        cos_theta_initial(1,i)=dot(H(i,:,1),gravity(1,:))/(sqrt(dot(H(i,:,1),H(i,:,1)))*sqrt(dot(gravity(1,:),gravity(1,:))));
        cos_theta(1,i)=dot(H(i,:,t),gravity(1,:))/(sqrt(dot(H(i,:,t),H(i,:,t)))*sqrt(dot(gravity(1,:),gravity(1,:))));
       if 0.34<sqrt(dot(H(i,:,t),H(i,:,t)))<0.68
          if sign(cos_theta(1,i))==sign(cos_theta_initial(1,i))
             if sqrt(dot((H_target(1,:)-Tar(i,:,t)),(H_target(1,:)-Tar(i,:,t))))<=delta
              miu=exp(c/((sqrt(dot(H(i,:,t),H(i,:,t)))-0.34)*(sqrt(dot(H(i,:,t),H(i,:,t)))-0.68)));
              X(i,:,t+1)=(round(r5)+(1-round(r5))*abs(normrnd(0,(1/exp(t)))))*X(i,:,t)+rand*lamda*miu*(rand*H_target(1,:)+rand*(X_best_global(t,:)-X(i,:,t)))+(1-round(r5))*0.0001*sign(randn)*round(rand);
                 for j=1:d
                     if X(i,j,t+1)>=ub
                       X(i,j,t+1)=ub;
                     else
                       X(i,j,t+1)=max(lb,X(i,j,t+1));
                     end
                 end
              fit(1,i)=cec17_func((X(i,:,t+1))',funclabel);
             else
              miu=exp(c/((sqrt(dot(H(i,:,t),H(i,:,t)))-0.34)*(sqrt(dot(H(i,:,t),H(i,:,t)))-0.68)));
              X(i,:,t+1)=(round(r5)+(1-round(r5))*abs(normrnd(0,(1/exp(t)))))*X(i,:,t)+rand*lamda*miu*(X_best_global(t,:)-X(i,:,t));
                for j=1:d
                     if X(i,j,t+1)>=ub
                       X(i,j,t+1)=ub;
                     else
                       X(i,j,t+1)=max(lb,X(i,j,t+1));
                     end
                 end

              fit(1,i)=cec17_func((X(i,:,t+1))',funclabel);
             end
          else
          X(i,:,t+1)=X(i,:,t);
                 for j=1:d
                     if X(i,j,t+1)>=ub
                       X(i,j,t+1)=ub;
                     else
                       X(i,j,t+1)=max(lb,X(i,j,t+1));
                     end
                 end
          fit(1,i)=cec17_func((X(i,:,t+1))',funclabel);
          end
       else
       sigma=power(((1.3293*sin(0.5*beta*pi))/(0.9064*beta*power(2,0.25))),1/beta);
       levy=0.03*((randn*sigma)/power(abs(randn),1/beta));
       R=unifrnd(-1,1,1,d);
       X(i,:,t+1)=(round(r5)+(1-round(r5))*abs(normrnd(0,(1/exp(t)))))*X(i,:,t)+rand*lamda*levy*R*(1/sqrt(dot(R,R)))+(1-round(r5))*0.0001*sign(randn)*round(rand);
                 for j=1:d
                    if X(i,j,t+1)>=ub
                       X(i,j,t+1)=ub;
                     else
                       X(i,j,t+1)=max(lb,X(i,j,t+1));
                    end
                 end
       fit(1,i)=cec17_func((X(i,:,t+1))',funclabel);
       end
        R_2(i,:)=unifrnd(-1,1,1,d);
       Tar(i,:,t+1)=(1/sqrt(dot(R_2(1,:),R_2(1,:))))*R_2(1,:);
       H(i,:,t+1)=(h_min+rand*(h_max-h_min))*Tar(i,:,t+1);
    end
    fit_rank(1,:)=sort(fit(1,:));
     for l=1:5
       for ii=1:n
            if fit(1,ii)==fit_rank(1,l);
            index_rank(1,l)=ii;
            end
       end
       X_best_rank(l,:)=X(index_rank(1,l),:,t);%五个最好的位置
     end
    for i=26:50
     u=unifrnd(-1,1,1,d);
     X(i,:,t)=(round(r5)+(1-round(r5))*abs(normrnd(0,(1/exp(t)))))*X_best_rank(randi(5),:)+lamda_find*randn*(u/sqrt(dot(u,u)))+(1-round(r5))*0.0001*sign(randn)*round(rand); 
     fit(1,i)=cec17_func((X(i,:,t+1))',funclabel);
    end
     [fit_best(1,t+1),index(1,t+1)]=min(fit(1,:));
    X_best(t+1,:)=X(index(1,t+1),:,t+1);
      if fit_best(1,t+1)<fit_best_global(1,t)
       fit_best_global(1,t+1)=fit_best(1,t+1);
       X_best_global(t+1,:)=X_best(t+1,:);
      else
        fit_best_global(1,t+1)=fit_best_global(1,t);
        X_best_global(t+1,:)=X_best_global(t,:);
      end
end

