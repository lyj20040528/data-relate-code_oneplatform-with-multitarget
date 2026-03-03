%CATEG_RND  Draws samples from a given one dimensional discrete distribution
%
% Syntax:
%   C = CATEG_RND(P,N)
%
% Author:
%   Simo S鋜kk? 2002
%
% In:
%   P - Discrete distribution, which can be a numeric array
%       of probabilities or a cell array of particle structures,
%       whose weights represent the distribution.
%   N - Number of samples (optional, default 1)
%
% Out:
%   C - Samples in a Nx1 vector
%
% Description:
%   Draw random category

% History:
%    3.12.2002  The first version
%
% Copyright (C) 2002 Simo S鋜kk?%               2008 Jouni Hartikainen
%
% $Id: categ_rnd.m,v 1.1.1.1 2003/09/15 10:54:35 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function C = categ_rnd(P,N)%%选取得个数
  %该函数完成按概率随机抽取的工作
  % Check arguments
  %
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    N = [];
  end
 
  if isempty(N)
    N = 1;
  end
  %
  % If particle structures are given
  %
  if iscell(P) & isfield([P{:}],'W')%iscell(p)：判断p是不是cell   数组。isfield([P{:}],'W'判断结构体里面有没有字段W
      tmp = [P{:}];
      P = [tmp.W]; %也就是如果P包含了除了目标出现概率以外的信息，那么做信息过滤，只要概率这一块
  end
  %
  % Draw the categories
  %
  C = zeros(1,N);
  P = P(:) ./ sum(P(:));%归一化 举例p=[1 1 2 1]，归一化后p=[0.2 0.2 0.4 0.2]
  P = cumsum(P);%进行累积求和，所以向量P最后的元素就是1,p=[0.2 0.4 0.8 1]
  for i=1:N
    C(i) = min(find(P > rand));%find(P > rand)找到比随机生成的[0 1)中的数大的P的最小元素P{i}。
  end%其实结合上面例子你会发现，返回的C(i)就是i，其实也就是按照出现概率进行抽样。
