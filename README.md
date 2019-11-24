# kron_generate
对Graph500的kron图生成器部分进行了复写，可以生成任意点数的32位二进制文件

本程序是根据graph500的octave程序简化后改编而来，源程序是：

function ijw = kronecker_generator (SCALE, edgefactor)
%% Generate an edgelist according to the Graph500 parameters.  In this
%% sample, the edge list is returned in an array with three rows,
%% where StartVertex is first row, EndVertex is the second row, and
%% Weight is the third row.  The vertex labels start at zero.
%%  
%% Example, creating a sparse matrix for viewing:
%%   ijw = kronecker_generator (10, 16);
%%   G = sparse (ijw(1,:)+1, ijw(2,:)+1, ones (1, size (ijw, 2)));
%%   spy (G);
%% The spy plot should appear fairly dense. Any locality
%% is removed by the final permutations.

  %% Set number of vertices.
  N = 2^SCALE;

  %% Set number of edges.
  M = edgefactor * N;

  %% Set initiator probabilities.
  [A, B, C] = deal (0.57, 0.19, 0.19);

  %% Create index arrays.
  ijw = ones (3, M);
  %% Loop over each order of bit.
  ab = A + B;
  c_norm = C/(1 - (A + B));
  a_norm = A/(A + B);

  for ib = 1:SCALE,
    %% Compare with probabilities and set bits of indices.
    ii_bit = rand (1, M) > ab;
    jj_bit = rand (1, M) > ( c_norm * ii_bit + a_norm * not (ii_bit) );
    ijw(1:2,:) = ijw(1:2,:) + 2^(ib-1) * [ii_bit; jj_bit];
  end

  %% Generate weights
  ijw(3,:) = unifrnd(0, 1, M);

  %% Permute vertex labels
  p = randperm (N);
  ijw(1:2,:) = p(ijw(1:2,:));

  %% Permute the edge list
  p = randperm (M);
  ijw = ijw(:, p);

  %% Adjust to zero-based labels.
  ijw(1:2,:) = ijw(1:2,:) - 1;

endfunction

由于只计算三角形时不需要第三个参数，原文件可以改成更简单的m程序

function ijw = kronecker_generator (SCALE, edgefactor)

  %% Set number of vertices.
  N = 2^SCALE;

  %% Set number of edges.
  M = edgefactor * N;

  %% Set initiator probabilities.
  [A, B, C] = deal (0.57, 0.19, 0.19);

  %% Create index arrays.
  ijw = ones (2, M);
  %% Loop over each order of bit.
  ab = A + B;
  c_norm = C/(1 - (A + B));
  a_norm = A/(A + B);

  for ib = 1:SCALE,
    %% Compare with probabilities and set bits of indices.
    ii_bit = rand (1, M) > ab;
    jj_bit = rand (1, M) > ( c_norm * ii_bit + a_norm * not (ii_bit) );
    ijw = ijw + 2^(ib-1) * [ii_bit; jj_bit];
  end

  %% Permute vertex labels
  p = randperm (N);
  ijw = p(ijw);

  %% Permute the edge list
  p = randperm (M);
  ijw = ijw(:, p);

  %% Adjust to zero-based labels.
  ijw = ijw - 1;
