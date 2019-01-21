clear all; close all; clc;

% Test 1
mat = [[0, 1]; [1, 0]];
pertIdx = 1;
obsIdx = [2, 1];

truePertOrders = [1, 2];
test1 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 2
mat = [[0, 1]; [1, 0]];
pertIdx = 2;
obsIdx = [1, 2];

truePertOrders = [2, 1];
test2 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 3
mat = [[0, 0]; [1, 0]];
pertIdx = [2, 1];
obsIdx = [1, 2];

truePertOrders = [[Inf, 1]; [1, 2]];
test3 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 4
mat = [[0, 1]; [1, 0]];
pertIdx = [1, 2];
obsIdx = 2;

truePertOrders = [[1, 2]; [NaN, 1]];
test4 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)

% Test 5
mat = [[0, 1]; [1, 0]];
pertIdx = [2, 1];
obsIdx = [];

truePertOrders = [[NaN, 1]; [1, NaN]];
test5 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 6
mat = [[0, 0, 1]; [1, 0, 1]; [0, 1, 0]];
pertIdx = [3, 1, 2];
obsIdx = [1, 3, 2];

truePertOrders = [[2, 2, 1]; [1, 2, 3]; [3, 1, 2]];
test6 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 7
mat = [[0, 1, 1, 0]; [1, 0, 1, 1]; [1, 1, 0, 0]; [0, 1, 1, 0]];
pertIdx = [3, 1, 4, 2];
obsIdx = [4, 3, 1, 2];

truePertOrders = [[2, 2, 1, 2]; [1, 2, 2, 3]; [3, 2, 3, 1]; [2, 1, 2, 2]];
test7 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 8
mat = [[0, 1, 1, 0]; [1, 0, 1, 1]; [1, 1, 0, 0]; [0, 1, 1, 0]];
pertIdx = [3, 1, 4, 2];
obsIdx = [1, 2];

truePertOrders = [[2, 2, 1, NaN]; [1, 2, NaN, NaN]; [3, 2, NaN, 1]; [2, 1, NaN, NaN]];
test7 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)


% Test 9
mat = [[0, 1, 1, 0]; [1, 0, 1, 1]; [1, 1, 0, 0]; [0, 1, 1, 0]];
pertIdx = [3, 1, 4, 2];
obsIdx = [3, 4];

truePertOrders = [[NaN, NaN, 1, 2]; [1, NaN, 2, 3]; [NaN, NaN, 2, 1]; [NaN, 1, 2, 2]];
test7 = nnz(TruePertOrders(mat, pertIdx, obsIdx) ~= truePertOrders)
