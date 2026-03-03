%GAUSS_PDF  Multivariate Gaussian PDF
%
% Syntax:
%   P = GAUSS_PDF(X,M,S)
%
% Author:
%   Simo S鋜kk? 2002
%
% In:
%   X - Dx1 value or N values as DxN matrix
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%
% Out:
%   P - Probability of X. 
%   
% Description:
%   Calculate values of PDF (Probability Density
%   Function) of multivariate Gaussian distribution
%
%    N(X |燤, S)
%
%   Function returns probability of X in PDF. If multiple
%   X's or M's are given (as multiple columns), function
%   returns probabilities for each of them. X's and M's are
%   repeated to match each other, S must be the same for all.
%
% See also:
%   GAUSS_RND

% History:
%   14.05.2003  Returns also the energy
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo S鋜kk?
%
% $Id: gauss_pdf.m,v 1.1.1.1 2003/09/15 10:54:34 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [P,E] = gauss_pdf(X,M,S)
%右边 X:当前量测。M：IM。S：IS
%左边 P：
%IM = H*X;为预期观测值。 IS = (R + H*P*H')。
  if size(M,2) == 1
    DX = X-repmat(M,1,size(X,2));%在本代码环境里，接收到的量测和预测的量测肯定是同尺寸的，这句话做的是尺寸的统一。
    E = 0.5*sum(DX.*(S\DX),1);%S\DX也就是S的逆矩阵乘DX。.*是逐个元素相乘,对于一维量测，也就是得到0.5sum(IS,1)
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  elseif size(X,2) == 1
    DX = repmat(X,1,size(M,2))-M;  
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  else
    DX = X-M;  
    E = 0.5*DX'*(S\DX);%DX'为其转置
    d = size(M,1);%对应高斯分布的维数
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));%%det求矩阵的行列式，matlab中log默认e为底
    P = exp(-E);
  end
