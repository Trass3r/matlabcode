function rpc
if 0
	tic;
	M = importdata('D:\\ttt.txt');
	CKPs = importdata('D:\\tttchk.txt');
	toc
else
	fid = fopen('D:\\ttt.bin', 'rb');
	fid2 = fopen('D:\\tttchk.bin', 'rb');
	M = fread(fid, inf, '*double');
	M = reshape(M, 5, [])';
	CKPs = fread(fid2, inf, '*double');
	CKPs = reshape(CKPs, 5, [])';
	fclose(fid);
	fclose(fid2);
end
% line sample    X Y Z

numGCPs = size(M, 1);
numCKPs = size(CKPs, 1);
%numPtsPerLayer = size(A, 1);
numLayers = length(unique(M(:,5)));
%numLines = numPtsPerLayer / numSamples;

if 0
	uniqueLines = unique(M(:,1));
	uniqueSamples = unique(M(:,2));
	uniqueHeights = unique(M(:,5));

	%numSamples = length(0:100:20063);
	numSamples = length(uniqueSamples);

	uniqueLines = uniqueLines(1 : round(numLines/10) : numLines);
	uniqueSamples = uniqueSamples(1 : round(numSamples/10) : numSamples);
	uniqueHeights = uniqueHeights(1 : round(numLayers/5) : numLayers);
end
%pick = 
%A = M(pick, :);

%% normalize data
mins = min(min(M), min(CKPs));
maxs = max(max(M), max(CKPs))
(mins + maxs)/2
(mins - maxs)/2
M2 = bsxfun(@minus,   M,  (maxs+mins)/2);
M2 = bsxfun(@rdivide, M2, (maxs-mins)/2);
CKPs2 = bsxfun(@minus,   CKPs,  (maxs+mins)/2);
clear CKPs;
CKPs2 = bsxfun(@rdivide, CKPs2, (maxs-mins)/2);

clear mins maxs;
%M2 = normalize(M);
%CKPs2 = normalize(CKPs);

%%
if 0
	trainIdxs = [];
	chkIdxs = [];
	for i = 0:numLayers-1
		trainIdxs = [trainIdxs (1:2:numPtsPerLayer + i*numPtsPerLayer)];
		chkIdxs   = [chkIdxs   (2:2:numPtsPerLayer + i*numPtsPerLayer)];
	end
	assert(isempty(intersect(trainIdxs, chkIdxs)));
	assert(isequal(union(trainIdxs, chkIdxs), 1:size(M,1)));

	T = M(trainIdxs, :);
	C = M(chkIdxs, :);
	return;
end

%A = M(M(:,5) == 1000, :);
%A = M(mod(M(:, 1), 200) == 0 & mod(M(:, 2), 200) == 0, :);
%A = A(1:10:size(A,1), :);
A = sortrows(M2, [1 2 -5]); % sort by line and sample

if numGCPs < 4100
%figure;
scatter3(A(:,3), A(:,4), A(:,5), 3, 'r');
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
%hold on;
%scatter3(CKPs2(:,3), CKPs2(:,4), CKPs2(:,5), 3, 'b');
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
%hold off;

clear A;

%% display lines
meanterrainheight = 1.5;%2200;
% one line per column
% line(X,Y,Z)
% X = start0  start1
%     middle0 middle1
%     end0    end1
for i = 0:numGCPs/numLayers-1
	X(1, i+1) = -0.5 * A(i*numLayers+1, 2); % TODO: determine why I need to flip it here
	Y(1, i+1) = A(i*numLayers+1, 1);
	Z(1, i+1) = meanterrainheight;
	for j = 0:numLayers-1
		X(j+2, i+1) = A(i*numLayers+1+j, 3);
		Y(j+2, i+1) = A(i*numLayers+1+j, 4);
		Z(j+2, i+1) = A(i*numLayers+1+j, 5);
	end
end
line(X,Y,Z);
%line(X(2:end, :), Y(2:end, :), Z(2:end, :));
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
%set(gca, 'ZDir', 'reverse');
end

%% setup B matrices
Bs = zeros(numGCPs, 39);
Bl = zeros(numGCPs, 39);
X = M2(:, 3);
Y = M2(:, 4);
Z = M2(:, 5);
%                     1          2  3  4   5      6    7     8     9     10                                         15                        18
Bl(:, 1:20) = [ones(numGCPs, 1), X, Y, Z, X.*Y, X.*Z, Y.*Z, X.*X, Y.*Y, Z.*Z, X.*Y.*Z, X.*Y.*Y, X.*Z.*Z, Y.*X.*X, Y.*Z.*Z, Z.*X.*X, Z.*Y.*Y, X.^3, Y.^3, Z.^3];
Bs(:, 1:20) = Bl(:, 1:20);
Bl(:, 21:end) = bsxfun(@times, -M2(:, 1), Bl(:, 2:20));
Bs(:, 21:end) = bsxfun(@times, -M2(:, 2), Bs(:, 2:20));

numCoeffs = 20; % 7, 10, 20
B = Bs(:, 1:numCoeffs);

%% fit polynomial
polynomWhich = 2; % 1 = line 2 = sample
%cond(B' * B)
condB = cond(B)
%polynomCoeffs = (B' * B) \ (B' * M2(:, polynomWhich)); %  + 0.001*eye(numCoeffs)
polynomCoeffs = B \ M2(:, polynomWhich);
residuals = B * polynomCoeffs;
residuals = residuals - M2(:, polynomWhich);
RMSE = sqrt(sum(residuals .* residuals)/numGCPs);
polyerrors = [min(residuals), max(residuals), RMSE, RMSE * max(M(:,polynomWhich)) / 2] % error min, max, rmse, rmse in pixels

%% fit RPCs
BsInfo = [rank(Bs) cond(Bs) size(Bs)]
%cond(Bs' * Bs)
BlInfo = [rank(Bl) cond(Bl) size(Bl)]
%cond(Bl' * Bl)
%linesCoeffs   = (Bl' * Bl) \ (Bl' * M2(:, 1));
%samplesCoeffs = (Bs' * Bs) \ (Bs' * M2(:, 2));
linesCoeffs   = Bl \ M2(:, 1); % uses QR
samplesCoeffs = Bs \ M2(:, 2);
%[tmp,idxs] = sort(abs(linesCoeffs));
%linesCoeffs(idxs(1:10)) = 0;
%[tmp,idxs] = sort(abs(samplesCoeffs));
%samplesCoeffs(idxs(1:10)) = 0;

%% Tikhonov
range = 0:0.002:1;
% 0.002 seems best
if 0
	i = 1;
	for h=range
		linesCoeffs2  = [Bl; sqrt(h)*eye(39)] \ [M2(:, 1); zeros(39,1)]; % uses QR
		samplesCoeffs2= [Bs; sqrt(h)*eye(39)] \ [M2(:, 2); zeros(39,1)];
		
		residualsL = Bl * linesCoeffs2;
		residualsS = Bs * samplesCoeffs2;

		residualsL = M2(:, 1) - residualsL;
		residualsS = M2(:, 2) - residualsS;

		RMSEl = sqrt(sum(residualsL .* residualsL)/numGCPs);
		RMSEs = sqrt(sum(residualsS .* residualsS)/numGCPs);
		minsL(i) = min(abs(residualsL));
		maxsL(i) = max(abs(residualsL));
		RMSEsL(i) = RMSEl;
		RMSEs2L(i) = RMSEl * max(M(:,1)) / 2;
	
		minsS(i) = min(abs(residualsS));
		maxsS(i) = max(abs(residualsS));
		RMSEsS(i) = RMSEs;
		RMSEs2S(i) = RMSEs * max(M(:,2)) / 2;
		i = i+1
	end
	figure;
	plot(range, [minsL; maxsL]);
end


%% svd variant
if 0
	% S contains the singular values
	[U,S,V] = svd(Bl, 'econ');
	% inverse of diagonal
	for i = 1:size(S,1)
		S(i,i) = 1 / S(i,i);
	end
	%S(end, end) = 0;
	linesCoeffs = V * S * U' * M2(:, 1);
end

%% now RPC iterations
%residualsS = zeros(numGCPs, 1);
%residualsL = zeros(numGCPs, 1);
%for i = 1:numCKPs
%	X = M2(i, 3);
%	Y = M2(i, 4);
%	Z = M2(i, 5);
%	BB = [1, X, Y, Z, X*Y, X*Z, Y*Z, X*X, Y*Y, Z*Z, X*Y*Z, X*Y*Y, X*Z*Z, Y*X*X, Y*Z*Z, Z*X*X, Z*Y*Y, X^3, Y^3, Z^3];
%	BBl = [BB, -M2(i, 1) * BB(2:end)];
%	BBs = [BB, -M2(i, 2) * BB(2:end)];

	% first we just save the result X*beta
%	residualsL(i) = Bl * linesCoeffs;
%	residualsS(i) = Bs * samplesCoeffs;

%	BB * samplesCoeffs
%	BB * samplesCoeffs - CKPs2(i, 2)
%	(BB * samplesCoeffs - CKPs2(i, 2)) * 20100 / 2
%end
%pack; % try to defragment
delta = 1;
while delta > 0.0001
	
	residualsL = Bl * linesCoeffs;
	residualsS = Bs * samplesCoeffs;

	residualsL = M2(:, 1) - residualsL;
	residualsS = M2(:, 2) - residualsS;

	RMSEl = sqrt(sum(residualsL .* residualsL)/numGCPs);
	RMSEs = sqrt(sum(residualsS .* residualsS)/numGCPs);
	linErrorsL = [min(abs(residualsL)), max(abs(residualsL)), RMSEl, RMSEl * max(M(:,1)) / 2] % error min, max, rmse, rmse in pixels
	linErrorsS = [min(abs(residualsS)), max(abs(residualsS)), RMSEs, RMSEs * max(M(:,2)) / 2]

	% now the real ones
	residualsL = (Bl(:, 1:20) * linesCoeffs(1:20)) ./ (1 + Bl(:, 2:20) * linesCoeffs(21:end));
	residualsS = (Bs(:, 1:20) * samplesCoeffs(1:20)) ./ (1 + Bs(:, 2:20) * samplesCoeffs(21:end));

	residualsL = M2(:, 1) - residualsL;
	residualsS = M2(:, 2) - residualsS;

	RMSEl = sqrt(sum(residualsL .* residualsL)/numGCPs);
	RMSEs = sqrt(sum(residualsS .* residualsS)/numGCPs);
	errorsL = [min(abs(residualsL)), max(abs(residualsL)), RMSEl, RMSEl * max(M(:,1)) / 2] % error min, max, rmse, rmse in pixels
	errorsS = [min(abs(residualsS)), max(abs(residualsS)), RMSEs, RMSEs * max(M(:,2)) / 2]

	
%	Wl = diag([1 1./(Bl(:, 2:20) * linesCoeffs(21:end))]);
%	Ws = diag([1 1./(Bs(:, 2:20) * samplesCoeffs(21:end))]);

	Wl = 1 + Bl(:, 2:20) * linesCoeffs(21:end); % column vector
	Ws = 1 + Bs(:, 2:20) * samplesCoeffs(21:end);
	plot(1:numGCPs, Wl, 1:numGCPs, Ws); legend Wl Ws; xlabel GCPs; ylabel denominator;
	% draw line for 0 crossing and all the poles
	inds = find(abs(residualsL) > 1.1);
	line([0 inds'; numGCPs inds'], [0 3.5*ones(1, length(inds)); 0 -1.5*ones(1, length(inds))], 'LineWidth', 2);

	drawnow;
	% for reference: theoretical matrix based approach
%	samplesCoeffs  = (Ws*Bs) \ (Ws*residualsS);

% TODO: modify Bl directly instead of here
% this is weight matrix "imported" into B
%	delta = bsxfun(@rdivide, Bl, Wl) \ (residualsL ./ Wl); % uses QR
%	linesCoeffs = linesCoeffs + delta;
%	delta = bsxfun(@rdivide, Bs, Ws) \ (residualsS ./ Ws);
%	samplesCoeffs = samplesCoeffs + delta;

	[Q1L, R1L] = qr(bsxfun(@rdivide, Bl, Wl), 0); % economic QR
	[Q1S, R1S] = qr(bsxfun(@rdivide, Bs, Ws), 0); % economic QR
	linesCoeffs   = (R1L + inv(R1L')) \ (Q1L'* (M2(:, 1) ./ Wl) + inv(R1L')*linesCoeffs);
	samplesCoeffs = (R1S + inv(R1S')) \ (Q1S'* (M2(:, 2) ./ Ws) + inv(R1S')*samplesCoeffs);
	
%	delta = sqrt(sum(delta .* delta));
end

return;

%% check checkpoints
chkLines = zeros(numCKPs, 1);
for i = 1:numCKPs
	X = CKPs2(i, 3);
	Y = CKPs2(i, 4);
	Z = CKPs2(i, 5);
	BB = [1, X, Y, Z, X*Y, X*Z, Y*Z, X*X, Y*Y, Z*Z, X*Y*Z, X*Y*Y, X*Z*Z, Y*X*X, Y*Z*Z, Z*X*X, Z*Y*Y, X^3, Y^3, Z^3];
	
%	CKPs2(i, 1)
%	BB * linesCoeffs
%	(BB * linesCoeffs - CKPs2(i, 1)) * 15000 / 2

%	errors(i) = BB(1:numCoeffs) * linesCoeffs - CKPs2(i, 1);
	chkLines(i) = BB(1:numCoeffs) * polynomCoeffs;

%	BB * samplesCoeffs
%	BB * samplesCoeffs - CKPs2(i, 2)
%	(BB * samplesCoeffs - CKPs2(i, 2)) * 20100 / 2
end
errors = chkLines - CKPs2(:, polynomWhich);
sqrt(sum(errors .* errors)/numCKPs)
ans * max(M(:,polynomWhich)) / 2 % error in pixels

if numGCPs < 4100
figure;
scatter3(A(:,3), A(:,4), A(:,5), 60, A(:,polynomWhich), '.');
hold on;
scatter3(CKPs2(:,3), CKPs2(:,4), CKPs2(:,5), 60, chkLines, '.');
colorbar;
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
end

%%
return;
figure;
A = M(M(:,5) == 0, :);
B = M(M(:,5) == 1000, :);
A = A(1:100:size(A,1), :);
B = B(1:100:size(B,1), :);


meanterrainheight = 2200;
% one line per column
% line(X,Y,Z)
% X = start0  start1
%     middle0 middle1
%     end0    end1
for i = 0:numGCPs/7-1
	X(1, i+1) = A(i*7, 1);
	Y(1, i+1) = A(i*7, 2);
	Z(1, i+1) = meanterrainheight;
	for j = 0:6
		X(j+2, i+1) = A(i*7+j, 3);
		Y(j+2, i+1) = A(i*7+j, 4);
		Z(j+2, i+1) = A(i*7+j, 5);
	end
end
X = [A(:,1)'
	 A(:,3)'];
Y = [A(:,2)'
	 A(:,4)'];
Z = [zeros(1,size(A,1))
	(A(:,5)' + 100)];
line(X,Y,Z);
xlabel 'X'; ylabel 'Y'; zlabel 'Z';
set(gca, 'ZDir', 'reverse');


function M2 = normalize(M)
mins = min(M);
maxs = max(M);
(mins + maxs)/2
(mins - maxs)/2
M2 = bsxfun(@minus, M, mins);
M2 = bsxfun(@rdivide, M2, (maxs - mins) / 2);
%M2 = bsxfun(@rdivide, M2, maxs / 2);
M2 = M2 - 1;
%min(M2)
%max(M2)