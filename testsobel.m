function [ output_args ] = testsobel( input_args )
%TESTSOBEL Summary of this function goes here
%   Detailed explanation goes here


s = [1 2 1; 0 0 0; -1 -2 -1];


% These commands extract the horizontal edges from a raised pedestal.
A = zeros(10);
A(3:7,3:7) = ones(5);
subplot 121
surf(A);
subplot 122
H = conv2(A,s);
mesh(H);