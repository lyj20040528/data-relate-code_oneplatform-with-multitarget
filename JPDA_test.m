clear all
close all

sd = 0.05;            % 噪声标准差
dt = 0.01;            % 测量时间间隔
nsteps = 500;         % 总步数
ntargets = 2;         % 目标个数
T = (1:nsteps)*dt;    % 时间
states = cell(ntargets,nsteps);%储存真实目标轨迹和速度信息

F = [0 0 1 0;
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];

%x1 = [-1; -0.6;  1.0;  0.10];%这两个目标离得很近，所以很容易出现关联错误(对于NNDA和GNNDA)
%x2 = [-1; -0.3; 1.0; 0.02];%两轴坐标和两轴速度

x1 = [-1;0;1;0.1];%两轴坐标和两轴速度
x2 = [4;0;-1;0.1];%x1(t)=[-1+t ,0.1t ,1 ,0.1]
                  %x2(t)=[4-t ,0.1t ,-1 ,0.1]
%jiaodian = [1.5;0.25];%两起始点减去交点
a = [-2.5; 2.5];    % 第一个向量：两个向量都是起点减去交点得到的
b = [-0.25; -0.25]; % 第二个向量
dot_product = a(1)*b(1) + a(2)*b(2);%计算点积
norm_a = norm(a);  %norm计算向量长度，即x、y分量平方和的根号
norm_b = norm(b);  %
cos_theta = dot_product / (norm_a * norm_b);
cos_theta = max(min(cos_theta, 1), -1);  % 强制限制在合法范围，截断在[-1,1]避免数值误差
angle_rad = acos(cos_theta);    % 弧度值（0-π）
jiaodu = rad2deg(angle_rad);  % 转换为角度（0-180°）
                  
%x1 = [-1;0;1;0.1];%两轴坐标和两轴速度
%x2 = [-1;0.5;1;-0.1];%这个相对第二套轨迹只改了起始点，但是运动过程相对夹角更小，就出现了错误关联(对于NNDA和GNNDA)

for k = 1:nsteps
    x1 = expm(F*dt)*x1;
    x2 = expm(F*dt)*x2;
    states{1,k} = x1;
    states{2,k} = x2;
end%两目标坐标和速度

% 真实轨迹
YY1 = zeros(2,nsteps);%2×nsteps
YY2 = zeros(2,nsteps);
for k = 1:nsteps
    YY1(:,k) = states{1,k}(1:2);
    YY2(:,k) = states{2,k}(1:2);
end%两目标坐标(即轨迹)

% 保存量测
Y = zeros(2,2,nsteps);%2×2×nsteps，一个周期两个量测
%Y = cell(2,nsteps);
Y1 = zeros(2,nsteps);%这两个用于画量测
Y2 = zeros(2,nsteps);

xmin = -1;
xmax = 4;
ymin = -0.2;
ymax = 0.7;%均匀分布生成杂波
pclutter = 0.03;
for k = 1:nsteps
    if rand < pclutter
        Y(:,1,k) = [xmin+(xmax-xmin)*rand; ymin+(ymax-ymin)*rand];%杂波项
    else
    Y(:,1,k) = states{1,k}(1:2) + sd*randn(2,1);%普通带噪量测
    end
    if rand < pclutter
        Y(:,2,k) = [xmin+(xmax-xmin)*rand; ymin+(ymax-ymin)*rand];
    else
    Y(:,2,k) = states{2,k}(1:2) + sd*randn(2,1);
    end
    Y1(:,k) = Y(:,1,k);%这个用来画量测
    Y2(:,k) = Y(:,2,k);
end
%
qx = 0.1;
qy = 0.1;
[A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);
H = [1 0 0 0; 0 1 0 0];
R = diag([sd^2 sd^2]);

%
X1 = cell(1,nsteps+1); P1 = cell(1,nsteps+1);
X2 = cell(1,nsteps+1); P2 = cell(1,nsteps+1);%1到501是因为跟踪前有初值

X1{1} = [0;0;0;0];  
P1{1} = diag([10 10 10 10]);%目标1

%X2{1} = [3;0;0;0]; %如果采用这个初值，可以不加只根据最大后验概率的事件进行更新的模块
X2{1} = [0;0;0;0]; %如果采用现在的初值，则要加上根据轨迹距离选择：是否只根据最大后验概率的事件进行更新的模块
P2{1} = diag([10 10 10 10]);%目标2

%X1{1} = [Y(:,1,1); 0; 0];  
%P1{1} = diag([10 10 10 10]);%目标1

%X2{1} = [Y(:,2,1); 0; 0]; 
%P2{1} = diag([10 10 10 10]);%目标2

% RTS
MM1 = zeros(4,nsteps); PP1 = zeros(4,4,nsteps);
MM2 = zeros(4,nsteps); PP2 = zeros(4,4,nsteps);
PG = 0.997;
m = 2;%目标个数
r = chi2inv(PG, m);%卡方分布
for i = 2:nsteps+1
    k = i-1;

    Xx1 = A * X1{i-1};
    Pp1 = A * P1{i-1} * A'+Q;       % 1的预测
    Xx2 = A * X2{i-1};
    Pp2 = A * P2{i-1} * A'+Q;       % 2的预测
%%
%计算各轨迹门内量测
v1 = cell(1,2);
S1 = R + H*Pp1*H';
v2 = cell(1,2);
S2 = R + H*Pp2*H';
d1 = [];%d1和d2用于存在轨迹跟踪门内的量测对应距离
d2 = [];
re1 = [];%re1和re2用于存在轨迹跟踪门内的量测编号
re2 = [];
%c1 = 0;
%c2 = 0;
  for j = 1:size(Y,2)%分别计算距离并且给门限用于筛选  
      z = Y(:,j,k);

    v1{j} = z - H*Xx1;%轨迹1分别计算和两个量测的统计距离和对应数据
    S1 = R + H*Pp1*H';
    d = v1{j}' / S1 * v1{j};%计算马氏距离，用于判断预测值和量测有多不一致,马氏距离＞0
    if d<r
        d1 = [d1 d];
        re1 = [re1 j];
    end
    v2{j} = z - H*Xx2;%轨迹2分别计算和两个量测的统计距离和对应数据
    S2 = R + H*Pp2*H';
    d = v2{j}' / S2 * v2{j};
     if d<r
        d2 = [d2 d];
        re2 = [re2 j];
    end
  end
  %%
  %求似然L(t,j)
  Pd = 1 - pclutter;                 %检测概率？
  Area = (xmax-xmin)*(ymax-ymin);
  c_den = 1/Area;            %杂波空间密度。杂波出现在观测空间内时的分布，并且是符合均匀分布。杂波似然
  
  % 似然
  %求每个目标t对每个门内量测j的高斯似然值L(t,j)j是量测个数，t是目标个数
  L = zeros(ntargets,size(Y,2));        % L(t,j) 求似然，是标量
  for j = 1:size(Y,2)%只有在门内才求似然，在门外似然就是0用于后续筛选
      z = Y(:,j,k);
      if any(re1==j)
          % Gaussian likelihood N(v;0,S)
          %Sj = S1{j};
          %vj = v1{j};
          %IM1 = z-v1{j};
          L(1,j) =  gauss_pdf(z,z-v1{j},S1);%这里用的是之前贺巍的似然
          %L(1,j) = exp(-0.5*(v1{j}'/S1{j}*v1{j})) / sqrt((2*pi)^2 * det(S1{j}));
      end
      if any(re2==j)
          %Sj = S2{j};
          %vj = v2{j};
          %IM2 = z-v2{j};
          L(2,j) =  gauss_pdf(z,z-v2{j},S2);
          %L(2,j) = exp(-0.5*(v2{j}'/S2{j}*v2{j})) / sqrt((2*pi)^2 * det(S2{j}));
      end
  end
%%
%算后验
  % Enumerate feasible joint events: a1,a2 in {0,1..mMeas}
  events = [];
  w = [];
  for a1 = 0:size(Y,2)
      for a2 = 0:size(Y,2)
          % constraint: a1,a2 cannot use same nonzero measurement
          if a1~=0 && a1==a2
              continue%跳出当前这个循环剩余的代码，也就是保证一个联合事件中量测不能共用
          end
          % gating feasibility
          if a1~=0 && L(1,a1)==0%如果a1不为0并且似然为0那么就代表该量测在轨迹1跟踪门外
              continue
          end
          if a2~=0 && L(2,a2)==0
              continue
          end
          % weight = prior * likelihood * clutter for unused measurements
          ww = 1;
          if a1==0%没有关联量测
              ww = ww*(1-Pd);
          else
              ww = ww*Pd*L(1,a1);
          end
          if a2==0
              ww = ww*(1-Pd);
          else
              ww = ww*Pd*L(2,a2);
          end

          used = [a1 a2];%关联量测组
          used = used(find(used~=0));%和下面一行代码效果一样
          %used = used(used~=0);%筛选出used中非零元素，比如used=[0 2],used(used~=0)得到2;used=[1 2],used(used~=0)得到[1 2]
          %nUnused = size(Y,2)-length(unique(used));%还剩多少量测没有被目标解释。
          nUnused = size(Y,2)-length(used);
          %unique函数返回的向量只保留used中重复元素中的一个。这个环境不需要unique函数
          %used=[1 1],则unique(used)=1
          ww = ww * (c_den^nUnused);%ww为当前事件的未归一化权重。(c_den^nUnused)是事件补上的杂波似然，因为先前只计算了正常量测似然
          events = [events; a1 a2];
          w = [w; ww];
      end
  end

  %if isempty(w) || sum(w)==0
      %X1{i} = Xx1;  P1{i} = Pp1;
      %X2{i} = Xx2;  P2{i} = Pp2;
  %else
%%
%归一化
      pE = w / sum(w);   %归一化权重
      beta = zeros(ntargets,size(Y,2));%目标t关联到量测j的后验概率
      beta0 = zeros(ntargets,1);%目标t漏检后验概率
      for e = 1:size(events,1)
          a1 = events(e,1); a2 = events(e,2);
          if a1==0
              beta0(1)=beta0(1)+pE(e); %如果说事件e中a1为0，那么该事件就包含轨迹1漏检，所以要加到漏检概率中
          else
              beta(1,a1)=beta(1,a1)+pE(e);%如果说事件e中a1不为0，那么该事件就包含轨迹1关联量测a1，所以要加到关联量测a1的概率中
          end
          
          if a2==0
              beta0(2)=beta0(2)+pE(e); 
          else 
              beta(2,a2)=beta(2,a2)+pE(e);
          end
      end
%%
%如果轨迹距离过近则选则最大后验的事件进行更新，而不是加权更新
      if norm(H*Xx1 - H*Xx2) < 3*sd%norm计算统计距离，如果太小就改成最大后验更新
          [~, emax] = max(pE);%获取pE中最大值索引
          a1m = events(emax,1); a2m = events(emax,2);
          beta(:,:) = 0;
          beta0(:) = 0;
          if a1m==0
              beta0(1)=1; 
          else
              beta(1,a1m)=1; 
          end
          if a2m==0
              beta0(2)=1; 
          else
              beta(2,a2m)=1;
          end
      end%把后验概率清零，然后把最大权值事件对应项置 1
%fprintf('total=%g\n', beta0(1)+sum(beta(1,:)));
%%
%对两个目标进行加权更新
      %更新目标1
      %S1c = R + H*Pp1*H';
      K1 = Pp1*H'/S1;
      vbar1 = zeros(2,1);
      C1 = zeros(2,2);   %
      %X11 = zeros(4,4);%zeros(4,4);
      %X11 = beta0(1)*Xx1*(Xx1)';
      for j = 1:size(Y,2)
          %vbar1 = vbar1 + beta(1,j)*(K1*v1{j}+Xx1)*(K1*v1{j}+Xx1)';%求和βjt(k)vj(k)。量测为j，目标为t
          vbar1 = vbar1 + beta(1,j)*v1{j};%求和βjt(k)vj(k)。量测为j，目标为t
          C1 = C1 + beta(1,j)*(v1{j}*v1{j}');
          %X11 = X11 + beta(1,j)*(Xx1 + K1*v1{j})*(Xx1 + K1*v1{j})';
      end
      X1{i} = Xx1 + K1*vbar1;
      %P1{i} = Pp1 - (1-beta0(1))*K1*S1*K1' + X11 - X1{i}*(X1{i})';%xj(k|k)=x(k|k-1)+Kvj
      %上面的形式也可也，两者形式一样
      P1{i} = Pp1 - (1-beta0(1))*K1*S1*K1' + K1*(C1 - vbar1*vbar1')*K1';%xj(k|k)=x(k|k-1)+Kvj。本式子推导无问题

      %更新目标2
      %S2c = R + H*Pp2*H';
      K2 = Pp2*H'/S2;
      vbar2 = zeros(2,1);
      C2 = zeros(2,2);
      %X22 = zeros(4,4);%zeros(4,4);
      for j = 1:size(Y,2)
          vbar2 = vbar2 + beta(2,j)*v2{j};
          C2 = C2 + beta(2,j)*(v2{j}*v2{j}');
      end
      X2{i} = Xx2 + K2*vbar2;
      P2{i} = Pp2 - (1-beta0(2))*K2*S2*K2' + K2*(C2 - vbar2*vbar2')*K2';
  %end

%%
    %RTS
    MM1(:,k)   = X1{i};
    PP1(:,:,k) = P1{i};
    MM2(:,k)   = X2{i};
    PP2(:,:,k) = P2{i};
end


%RTS smooth
[SM1, SP1, D1] = rts_smooth(MM1, PP1, A, Q);
[SM2, SP2, D2] = rts_smooth(MM2, PP2, A, Q);
%%
% 画图
M1 = zeros(2,nsteps);
M2 = zeros(2,nsteps);
for k=1:nsteps%只取坐标
    M1(:,k) = MM1(1:2,k);
    M2(:,k) = MM2(1:2,k);
end

% 画图
figure(1); 

h0(1)=plot(Y1(1,:),Y1(2,:),'+');  hold on            
h0(2)=plot(Y2(1,:),Y2(2,:),'+'); 
%h0(3)=plot(Y3(1,:),Y3(2,:),'+'); 杂波量测
set(h0,'color',[0.5 0.5 0.5]); 

h1=plot(M1(1,:),M1(2,:),'k-');           % KF 跟踪1
set(h1,'linewidth',2); 

h2=plot(M2(1,:),M2(2,:),'g-');           % KF 跟踪2
set(h2,'linewidth',2); 

h7=plot(YY1(1,:),YY1(2,:),'--');         % 真实轨迹1
set(h7,'color','r');
set(h7,'linewidth',2); 

h8=plot(YY2(1,:),YY2(2,:),'--');         % 真实轨迹2
set(h8,'color','b');
set(h8,'linewidth',2); 

h9=plot(M1(1,1),M1(2,1),'ko');           % KF 跟踪1的起始点
set(h9,'color','m');
set(h1,'markersize',10); 

h22=plot(M2(1,1),M2(2,1),'ko');           % KF 跟踪2的起始点
set(h22,'color','m');
set(h2,'markersize',10); 

%h5=plot(SM1(1,:), SM1(2,:), 'b--');      % RTS smooth 1
%set(h5,'linewidth',2); hold on

%h6=plot(SM2(1,:), SM2(2,:), 'r--');      % RTS smooth 2
%set(h6,'linewidth',2); hold on

xlabel('X轴 NNDA','fontsize',14)
ylabel('Y轴','fontsize',14)
xlim([-1.5 4.5])

hh=legend([h0(1),h1,h2,h7,h8,h9], ...
          '量测', '跟踪轨迹1', '跟踪轨迹2', ...
          '真实轨迹1','真实轨迹2','起始点');%因为h0(1)和h0(2)都是量测，我标注的时候标一个就行了
      
%hh=legend([h0(1),h1,h2,h3,h4,h5,h6], ...
%          '量测', '跟踪轨迹1', '跟踪轨迹2', ...
%          '真实轨迹1','真实轨迹2','RTS Smooth1','RTS Smooth2');%因为h0(1)和h0(2)都是量测，我标注的时候标一个就行了
set(hh,'fontsize',8)
%title(sprintf('JPDA跟踪结果'), 'FontSize',14);
title(sprintf('JPDA跟踪结果，夹角为 %.3f', jiaodu), 'FontSize',14);
%hold off;

figure(2)
h5=plot(SM1(1,:), SM1(2,:), 'k-');      % RTS smooth 1
set(h5,'linewidth',2); hold on

h6=plot(SM2(1,:), SM2(2,:), 'g-');      % RTS smooth 2
set(h6,'linewidth',2); 

h7=plot(YY1(1,:),YY1(2,:),'--');         % 真实轨迹1
set(h7,'color','r');
set(h7,'linewidth',2); 

h8=plot(YY2(1,:),YY2(2,:),'--');         % 真实轨迹2
set(h8,'color','b');
set(h8,'linewidth',2); 

h9(1)=plot(Y1(1,:),Y1(2,:),'+');            
h9(2)=plot(Y2(1,:),Y2(2,:),'+'); 
%h0(3)=plot(Y3(1,:),Y3(2,:),'+'); 杂波量测
set(h9,'color',[0.5 0.5 0.5]); 

xlabel('X轴','fontsize',14)
ylabel('Y轴','fontsize',14)
title(sprintf('JPDA平滑结果'), 'FontSize',14);
xlim([-1.5 4.5])

hh=legend([h9(1),h5,h6,h7,h8], ...
         '量测','轨迹1平滑', '轨迹2平滑','真实轨迹1','真实轨迹2');
set(hh,'fontsize',8)
%%
%RMSE，图3是跟踪的，图4是平滑的
X1 = zeros(1,nsteps);
Y1 = zeros(1,nsteps);
n1 = zeros(1,nsteps);
  for i=1:nsteps
      X1(i) = (M1(1,i) - YY1(1,i))^2;
      Y1(i) = (M1(2,i) - YY1(2,i))^2;
      n1(i) = sqrt((X1(i) + Y1(i)));%/(YY1(2,i)^2+YY1(1,i)^2);如果改成这个就把sqrt去掉
  end
X2 = zeros(1,nsteps);
Y2 = zeros(1,nsteps);
n2 = zeros(1,nsteps);
  for i=1:nsteps
      X2(i) = (M2(1,i) - YY2(1,i))^2;
      Y2(i) = (M2(2,i) - YY2(2,i))^2;
      n2(i) = sqrt((X2(i) + Y2(i)));%/(YY2(2,i)^2+YY2(1,i)^2);
  end
  
figure(3)
h1=plot(T,n1,'-');
set(h1,'color','k');
set(h1,'linewidth',2); 
hold on;
h2=plot(T,n2,'-');
set(h2,'color','g');
set(h2,'linewidth',2); 

xlabel('时间T','fontsize',14)
ylabel('RMSE误差','fontsize',14)
xlim([0 5])

hh=legend([h1,h2], '轨迹1RMSE', '轨迹2RMSE');
set(hh,'fontsize',12)
title(sprintf('跟踪结果RMSE'), 'FontSize',14);

X3 = zeros(1,nsteps);
Y3 = zeros(1,nsteps);
n3 = zeros(1,nsteps);
  for i=1:nsteps
      X3(i) = (SM1(1,i) - YY1(1,i))^2;
      Y3(i) = (SM1(2,i) - YY1(2,i))^2;
      n3(i) = sqrt((X3(i) + Y3(i)));%/(YY1(2,i)^2+YY1(1,i)^2);
  end
X4 = zeros(1,nsteps);
Y4 = zeros(1,nsteps);
n4 = zeros(1,nsteps);
  for i=1:nsteps
      X4(i) = (SM2(1,i) - YY2(1,i))^2;
      Y4(i) = (SM2(2,i) - YY2(2,i))^2;
      n4(i) = sqrt((X4(i) + Y4(i)));%/(YY2(2,i)^2+YY2(1,i)^2);
  end

figure(4)
h1=plot(T,n3,'-');
set(h1,'color','k');
set(h1,'linewidth',2); 
hold on;
h2=plot(T,n4,'-');
set(h2,'color','g');
set(h2,'linewidth',2); 

xlabel('时间T','fontsize',14)
ylabel('RMSE误差','fontsize',14)
xlim([0 5])

hh=legend([h1,h2], '轨迹1RMSE', '轨迹2RMSE');
set(hh,'fontsize',12)
title(sprintf('平滑结果RMSE'), 'FontSize',14);

