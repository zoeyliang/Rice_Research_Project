function ip = spinterpm(d,outputLength,z,y,levelseq,purgedata)
% SPINTERPM   Multi-linear interpolation (maximum-norm-based grid)
%    IP = SPINTERPM(N,D,Z,Y)  Computes the interpolated
%    values at [Y1, ..., YN] over the sparse grid at level N. Y
%    may be a double array for vectorized processing, with each row
%    representing one point. The sparse grid data must be given  as
%    an array Z containing the node weights (computed with
%    SPVALS). Note that the sparse grid is always normalized to
%    the unit cube [0,1]^D, i.e., if the weights have been computed
%    for a different domain, the values Y have to be rescaled
%    accordingly. (Internal function)
	
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.2
% Date   : February 3, 2006

% Change log:
% V1.0   : September 24, 2003
%          Initial revision
% V1.1   : January 24, 2006
%          Changed data types to operate on uint arrays
% V1.2   : February 3, 2006
%          Added droptol processing.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

ninterp = uint32(size(y,1));
ip = zeros(outputLength,ninterp);

% Get the number of levels
nlevels = uint32(size(levelseq,1));

if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
index2 = zeros(d,1,'uint32'); 
repvec = ones(d,1,'uint32');
level = ones(d,1,'uint8');
repeat = zeros(d,1);
		
for kl = 1:nlevels
	npoints = uint32(1);
	lval = uint8(0);
	nlevelzero = 0;
	for k = 1:d
		lval = levelseq(kl,k);
		level(k) = lval;
		if lval == 0
			repvec(k) = 3;
			nlevelzero = nlevelzero + 1;
		else
			repvec(k) = 2^uint32(lval);
		end
		repeat(k) = 0;
		npoints = npoints * repvec(k);
		if k > 1
			repvec(k) = repvec(k) * repvec(k-1);
		end
	end

	% Skip subgrids with all surpluses below droptol.
	if purge
		if purgedata(kl) == 0
			index = index + npoints;
			continue;
		end
	end
	
	nrepeats = 2^nlevelzero-1;

	yt = 0;
	
	k = uint32(1);
	repstep = 0;
	while k <= ninterp
		temp = 1;
		l = uint16(1);
		while l <= d
			lval = level(l);
			yt = y(k,l);
			
			% Compute the scaling factor and the array position of
			% the weight
			if lval == 0
				if repeat(l) == 0
					index2(l) = 1;
					temp = temp * (1 - 2 * abs(yt-0.5));
				else
					if yt == 1
						index2(l) = 2;
					else
						xp = floor(yt * 2) * 2;
						if xp == 0
							temp = temp * 2 * (0.5 - yt);
						else
							temp = temp * 2 * (yt - 0.5);
						end
						index2(l) = xp;
					end
				end
			elseif yt == 1
				temp = 0;
				break;
			else
				scale = 2^double(lval);
				xp = floor(yt * scale);
				temp = temp * ...
							 (1 - 2 * scale * abs( yt - (xp+0.5)/scale));
				index2(l) = xp;
			end
			l = l + 1;
			if temp == 0
				break;
			end
		end
		
		% If the scaling factor is not Zero, add the computed value
		if temp > 0
			index3 = index + index2(1);
			for l = 2:d
				index3 = index3 + repvec(l-1)*index2(l);
			end
			ip(:,k) = ip(:,k) + temp*z(:,index3);
		end
		if repstep == nrepeats
			k = k + 1;
			if nrepeats > 0
				for l = 1:d
					repeat(l) = 0;
				end
				repstep = 0;
			end
		else
			for l = 1:d
				if level(l) == 0
					repeat(l) = repeat(l) + 1;
					if repeat(l) > 1
						repeat(l) = 0;
					else
						break;
					end
				end
			end	
			repstep = repstep + 1;
		end
	end
	index = index + npoints;
end

