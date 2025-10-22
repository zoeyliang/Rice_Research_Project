function ipder = dctdiffcheb(n,z,outputLength,y)
% DCTDIFFCHEB  Compute the derivative of a univariate 
%    polynomial given by hierarchical surpluses at the
%    Cehbyshev-Gauss-Lobatto nodes) using FFT/IFFT. 
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : June 25, 2006

% Change log:
% V1.0   : Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------
% z = pagetranspose(z);%num2cell(z', 2);

ninterp = size(z,2);
ipder = zeros(outputLength,length(y));


% Extend grid values
zero_temp = cell(1,ninterp);
zero_temp(:) = {zeros(1,outputLength)};
if n == 3
	zext = [z(2,:); zero_temp; z(1,:); zero_temp];
else
    zext = cell((n-1)*2,ninterp);
    zext(:) = {zeros(1,outputLength)};
    zext(2:2:end,:) = [flipud(z); z];
end
zext = cell2mat(zext);

% Perform DCT via FFT
chebc = real(fft(zext))/(double(n)-1);

chebc(n,:) = chebc(n,:)/2;
chebc(1,:) = chebc(1,:)/2;

% Differentiate (in place)
% for k = 1:ninterp
	c = chebc(n-1,:);
	chebc(n-1,:) = 2*double(n-1)*chebc(n,:);
	chebc(n,:)   = 0;
	l = double(n) - 2;
	while l >= 2
	  cprev = c;
		c = chebc(l,:);
		chebc(l,:) = chebc(l+2,:) + 2*l*cprev;
		l = l - 1;
	end
	chebc(1,:) = chebc(3,:)/2 + c;
% end

% Normalize to interval [0,1]
chebc = chebc * 2;

nn = size(chebc,1);
chebc = mat2cell(chebc,ones(1,nn),outputLength*ones(1,ninterp));



% Evaluate using the Clenshaw recurrence formula
if ninterp ~= length(y)
  % Univariate case is slightly different; there is only
	% a single column of surplus values in z
	ninterp = length(y);
% 	for k = 1:ninterp
% 		c  = cell(1,ninterp);
%         c(:) = {zeros(1,outputLength)};
% 		c2 = cell(1,ninterp);
%         c2(:) = {zeros(1,outputLength)};
        c = zeros(outputLength,ninterp);
        c2 = zeros(outputLength,ninterp);
		l  = n - 2;
		y1 = 2 * y - 1;
		y2 = 2*y1;
		tmp = 0;
		while l >= 1
			tmp = c;
%             cellfun(@times,A,num2cell(b'),'uni',false)
			c = bsxfun(@times,c,y2') - c2 + repmat(cell2mat(chebc(l+1,:))',1,ninterp);
			c2 = tmp;
			l = l - 1;
		end
		ipder = bsxfun(@times,c,y1') - c2 + repmat(cell2mat(chebc(1,:))',1,ninterp);

else
  % Process multiple columns from multivariate case
% 	for k = 1:ninterp
        c = zeros(outputLength,ninterp);
        c2 = zeros(outputLength,ninterp);
		l  = n - 2;
		y1 = 2 * y - 1;
		y2 = 2*y1;
		tmp = 0;
		while l >= 1
			tmp = c;
%             cellfun(@times,A,num2cell(b'),'uni',false)
			c = bsxfun(@times,c,y2') - c2 + cell2mat(reshape(chebc(l+1,:),ninterp,1))';
			c2 = tmp;
			l = l - 1;
		end
		ipder = bsxfun(@times,c,y1') - c2 + cell2mat(reshape(chebc(1,:),ninterp,1))';
% 	end
end
% ipder = cell2mat(reshape(ipder,[],1))';
