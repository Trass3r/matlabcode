function testdct()
%TESTDCT Summary of this function goes here
%   Detailed explanation goes here

A = [
	156 144 125 109 102 106 114 121
	151 138 120 104  97 100 109 116
	141 129 110  94  87  91  99 106
	128 116  97  82  75  78  86  93
	114 102  84  68  61  64  73  80
	102  89  71  55  48  51  60  67
	 92  80  61  45  38  42  50  57
	 86  74  56  40  33  36  45  52
	];

plot(1:8, A(1,:), 1:8, dct1(A(1,:)'));

function A = dct1(A)
	N = length(A);
	I = bsxfun(@times, 0.5:1:N, (0:N-1)');
	M = cos(pi/N * I);
	A = M * A;
end

dct2(A)

function A = dct2(A)
	[M,N] = size(A);
	I = bsxfun(@times, 0.5:1:N, (0:N-1)');
	M = cos(pi/N * I);
	first = (M * A')';
	A = M * first;
end
end