%RTS_SMOOTH  Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,S] = RTS_SMOOTH(M,P,A,Q)
%
% In:
%   M - NxK matrix of K mean estimates from Kalman filter
%   P - NxNxK matrix of K state covariances from Kalman Filter
%   A - NxN state transition matrix or NxNxK matrix of K state
%       transition matrices for each step.
%   Q - NxN process noise covariance matrix or NxNxK matrix
%       of K state process noise covariance matrices for each step.
%
% Out:
%   M - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   D - Smoother gain sequence
%   
% Description:
%   Rauch-Tung-Striebel smoother algorithm. Calculate "smoothed"
%   sequence from given Kalman filter output sequence
%   by conditioning all steps to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   MM = zeros(size(m,1),size(Y,2));
%   PP = zeros(size(m,1),size(m,1),size(Y,2));
%   for k=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,k),H,R);
%     MM(:,k) = m;
%     PP(:,:,k) = P;
%   end
%   [SM,SP] = rts_smooth(MM,PP,A,Q);
%
% See also:
%   KF_PREDICT, KF_UPDATE

% Copyright (C) 2003-2006 Simo S鋜kk?%
% $Id: rts_smooth.m 109 2007-09-04 08:32:58Z jmjharti $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

%





function [M,P,D] = rts_smooth(M,P,A,Q)%右边的M是平滑前均值；右边的P是平滑前的状态协方差；A是状态转移矩阵；
%Q是过程噪声协方差，并且是离散化之后的。
%式子左边分别是：M和P分别是平滑后的状态均值和协方差；D是平滑增益



  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end%参数不够就进行不了RTS平滑

  %
  % Extend A and Q if they are NxN matrices
  %
  
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);%repmat函数，复制函数。repmat(A,[a,b,c])，将数组A在1维复制a次，2维复制b次，3维复制c次。
  end                             %A是状态转移矩阵，如果它时不变，那么A在时间维度就不会变化，所以此时size(A,3)==1。
  if size(Q,3)==1                 %此时在时间维做复制，是为了后续for循环调用
    Q = repmat(Q,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));%初始化平滑增益
  %这里的M和KLM_2D的M不一样。KLM_2D的M只取坐标，为了方便画图。卡尔曼平滑里面的M，包含了完整信息，即x、y轴的位置和速度
  
  for k=(size(M,2)-1):-1:1%-1是for循环每次循环的步长，意思是从size(M,2)-1开始，没循环一次减1，一直循环到k=1
    P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);%对应论文k+1时刻预测的协方差Pk+1
    D(:,:,k) = P(:,:,k) * A(:,:,k)' / P_pred;
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - A(:,:,k) * M(:,k));%均值平滑
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_pred) * D(:,:,k)';%协方差平滑
  end

