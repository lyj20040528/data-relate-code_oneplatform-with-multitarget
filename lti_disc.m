%LTI_DISC  Discretize LTI ODE with Gaussian Noise
%
% Syntax:
%   [A,Q] = lti_disc(F,L,Qc,dt)
%
% In:
%   F  - NxN Feedback matrix
%   L  - NxL Noise effect matrix        (optional, default identity)
%   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
%   dt - Time Step                      (optional, default 1)
%
% Out:
%   A - Transition matrix
%   Q - Discrete Process Covariance
%
% Description:
%   Discretize LTI ODE with Gaussian Noise. The original
%   ODE model is in form
%
%     dx/dt = F x + L w,  w ~ N(0,Qc)
%
%   Result of discretization is the model
%
%     x[k] = A x[k-1] + q, q ~ N(0,Q)
%
%   Which can be used for integrating the model
%   exactly over time steps, which are multiples
%   of dt.

% History:
%   11.01.2003  Covariance propagation by matrix fractions
%   20.11.2002  The first official version.
%
% Copyright (C) 2002, 2003 Simo S鋜kk?
%
% $Id: lti_disc.m 111 2007-09-04 12:09:23Z ssarkka $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [A,Q] = lti_disc(F,L,Q,dt)%定义函数lti_disc
%L:反映噪声如何进入系统。
  
  % Check number of arguments
  %nargin函数负责输出调用函数lti时，上面的四个参数给了几个。
  %eg:lti_disc(F)--nargin=1     eg:lti_disc(F,L)--nargin=2
  
  %nargin 仅用于判断函数被调用时传入的参数个数，MATLAB 按参数位置进行传递，
  %因此不能跳过中间参数，只能通过传入空数组并在函数内部进行判断来实现可选参数机制。
  
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    L = [];%matlab里面，我想要给一个参数a设为空，那么就是a=[]
  end
  if nargin < 3
    Q = [];
  end
  if nargin < 4
    dt = [];
  end%不用 else if，而是 一连串独立的if目的就是：当参数给得不够时，后面的参数能顺序地被补上空值，无参数输入除外。

  if isempty(L)
    L = eye(size(F,1));
  end
  if isempty(Q)
    Q = zeros(size(F,1),size(F,1));
  end
  if isempty(dt)
    dt = 1;
  end%这一段就是给空值的参数赋值，nargin代码段负责给没传的参数赋空值，isempty负责给空值的部分赋值

  %
  % Closed form integration of transition matrix
  %
  A = expm(F*dt);

  %
  % Closed form integration of covariance
  % by matrix fraction decomposition
  %
  n   = size(F,1);%matlab中，矩阵∈数组。size有以下两种用法：
                  %对于数组A，一个是size(A)，得到n维数组A的各维的维数组成的行向量，一个是size(A,k)得到数组A第k维的维数
                  %eg A为3×4的矩阵，size(A)=[3,4],size(A,1)=3
  Phi = [F L*Q*L'; zeros(n,n) -F'];%矩阵拼接。[a,b]或者[a b]：横向拼接。[;]：纵向拼接
  %F是连续时间状态转移矩阵，LQL'反映连续过程噪声在状态转移时的不确定性强度，会用于后续计算离散过程噪声协方差
  AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];%AB 是一个中间量，用于把连续时间过程噪声在连续状态中产生的影响，
                                         %转换为离散时间下作用于离散状态的过程噪声协方差的计算形式。
                                         %AB 是连续白噪声在 dt 内积分后，对离散状态的等效影响。
  Q   = AB(1:n,:)/AB((n+1):(2*n),:);
  
