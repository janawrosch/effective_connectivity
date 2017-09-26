function varargout = gplotcolor(A, xy, cmap, pct)
%GPLOTCOLOR Plot graph theory-style graph with colored lines
%
% gplotcolor(A, xy, cmap)
% gplotcolor(A, xy, cmap, pct)
% h = gplotcolor(A, xy, cmap)
%
% This function plots a directed graph.  It differs from gplot in two ways:
% - it indicates direction of edge through a counterclockwise curve (stole
%   this idea from gplotdc.m, FEX 19580)
% - it plots each edge as a separate line object, allowing differing
%   colors/styles/etc for each edge
%
% Input variables:
%
%   A:      n x n adjacency matrix 
%
%   xy:     n x 2 position coordinates for nodes
%
%   cmap:   nnz x 3 colormap array, where nnz is the number of nonzero
%           values of A.  Colors will correspond to edges described by
%           find(A), respectively. 
%
%   pct:    curvature fraction [0.04]
%
% Output variables:
%
%   h:      handles of lines used to plot each edge

% Copyright 2010 Kelly Kearney

if ~isequal(length(find(A)), size(cmap,1))
    error('Colors array should have same number of points as non-zero elements');
end

[i,j] = find(A);

if nargin < 4
    pct = 0.04;
end
% [ignore, p] = sort(max(i,j));
% i = i(p);
% j = j(p);

% cmap = cmap(p,:);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';

[X,Y] = makeCurved(X,Y,pct);

h = plot(X,Y);
set(h, {'color'}, num2cell(cmap,2));

if nargout == 1
    varargout{1} = h;
end


function [xc,yc] = makeCurved(x,y,pct)

nx = size(x,2);
[xc,yc] = deal(zeros(17,nx));

for ix = 1:nx
    xtmp = x(:,ix);
    ytmp = y(:,ix);
    pctn = pct;
    for k = 1:4
        n = length(xtmp);
        m = 2*n - 1;
        xtmp2 = zeros(1,m);
        ytmp2 = zeros(1,m);
        xtmp2(1:2:m) = xtmp;
        ytmp2(1:2:m) = ytmp;
        xtmp2(2:2:m-1) = 0.5*(xtmp(1:n-1)+xtmp(2:n))+pctn*diff(ytmp);
        ytmp2(2:2:m-1) = 0.5*(ytmp(1:n-1)+ytmp(2:n))-pctn*diff(xtmp);
        xtmp = xtmp2;
        ytmp = ytmp2;
        pctn = pctn * 0.5;
    end
    xc(:,ix) = xtmp;
    yc(:,ix) = ytmp;
    
end
    

% 
% N = length(x);
% if N < 2
%     X = x;
%     Y = y;
% else
%     M = 2*N-1;
%     X = zeros(1,M);
%     Y = zeros(1,M);
%     X(1:2:M) = x;
%     Y(1:2:M) = y;
%     X(2:2:M-1) = 0.5*(x(1:N-1)+x(2:N))+pct*diff(y);
%     Y(2:2:M-1) = 0.5*(y(1:N-1)+y(2:N))-pct*diff(x);
% end
% PCT = 0.5*pct;