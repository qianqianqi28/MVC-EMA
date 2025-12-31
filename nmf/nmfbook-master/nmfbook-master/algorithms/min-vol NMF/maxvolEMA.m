% This codes modified the code "minvolNMF.m" on the website https://gitlab.com/ngillis/nmfbook/-/tree/master/algorithms/min-vol%20NMF?ref_type=heads
% On the website, the code is for minimum volume constrained NMF.
% Here the code is for maximum volume constrained EMA.
% 
% min_{W,H} 1/2||X-WH||_F^2 + lambda/2 * det(W^TW) 
% 
% where W >= 0, H >= 0, 
% and sum-to-one constraints on W or H: 
% H^T e <= e (model 1), or 
% H e    = e (model 2), or 
% W^T e  = e (model 3, default), or 
% H^T e  = e (model 4). 
% 
% This is solved by optimizing alternatively over W and H, 
% using a projected fast gradient method; see 
% V. Leplat, A.M.S. Ang, N. Gillis, "Minimum-Volume Rank-Deficient 
% Nonnegative Matrix Factorizations", ICASSP 2019, May 12-17, 2019, 
% Brighton, UK. 
% Zhou, G., Xie, S., Yang, Z., Yang, J.-M., & He, Z. (2011). Minimum-volume-constrained
% nonnegative matrix factorization: Enhanced ability of learning parts. IEEE Transac-
% tions on Neural Networks, 22(10), 1626-1637.
% 
% ****** Input ****** 
% X : m-by-n matrix to factorize 
% r : factorization rank r 
% options: 
%   .model: = 1: H^T e  <=  e, cols of H sum to at most one (not wlog).  
%           = 2:   H e   =  e, rows of H sum to one (wlog). 
%           = 3:   W^T e =  e, columns of W sum to one (wlog) --> default 
%           = 4: H^T e   =  e, cols of H sum to one (not wlog) --> *new*  
%                There is a new numerical experiment using this last model;
%                see the folder "examples by chapter/Chapter 4 - Identifiability/minvolNMF_Moffet.m"
%                and also https://www.dropbox.com/s/iqp2saalnolfm2r/errataandmore_NMFbook.pdf?dl=0 
%   .lambda' will be used to set up lambda (default: 0.1), 
%       lambda = lambda' * ||X-WH||_F^2 / det( W^TW) 
%       where (W,H) is the initialization. 
%       This allows to balance the two terms in the objective (and also
%       make the algorithm unsensitive to scaling of the input matrix). 
%   .maxiter (default: 100) 
%   .target: it allows to define a target for the relative error
%       ||X-WH||_F/||X||_F: lambda will be automatically tuned to attemp
%       achieving this target value (default: no target value). It is a 
%       simple heuristic: if ||X-WH||_F/||X||_F > target, lambda is 
%       decreased, and if ||X-WH||_F/||X||_F < target, lambda is increased. 
%   .(W,H): initialization (default: use of SNPA) 
%   .display: =1 displays the iteration count, = 0 no display
%               (default: 0).  
%
% ****** Outut ****** 
% (W,H) : low-rank approximation W*H of X, with H>=0, W>=0 of volume,  
%         det(W^TW) is big if lambda < 0, and  
%         if model = 1: H^T e  <=  e, or 
%                  = 2:   H e   =  e, or 
%                  = 3:   W^T e =  e (default), 
%                  = 4:   H^T e =  e. 
% e     : evolution of the error 
%         1/2||X-WH||_F^2 + lambda/2 * det( W^TW) 
% err1  = evolution of ||X-WH||_F^2
% err2  = evolution of det( W^TW) 

function [W,H,e,err1,err2, lambda] = maxvolEMA(X,r,options);

if nargin <= 2
    options = [];
end
if ~isfield(options,'model')
    options.model=3;
end
if ~isfield(options,'lambda')
    options.lambda=0.1;
end
if ~isfield(options,'maxiter')
    options.maxiter=100;
end
if ~isfield(options,'inneriter')
    options.inneriter=10;
end
if isfield(options,'W') && isfield(options,'H') 
    W = options.W;
    H = options.H;
else
    if options.model == 1
        options.proj = 1;
    end
    [K,H] = SNPA(X,r,options); 
    W = X(:,K);
    optionsNNLS.init = H; 
    H = NNLS(W,X,optionsNNLS); 
    if length(K) < r
        warning('SNPA recovered less than r basis vectors.');
        warning('The data poins have less than r vertices.');
        r = length(K);
        fprintf('The new value of r is %2.0d.\n',r);
    end
end
if ~isfield(options,'display')
    options.display=1;
end 
% Normalization 
[W,H] = normalizeWH(W,H,options.model,X); 
% Initializations
normX2 = sum(X(:).^2);
normX = sqrt(normX2); 
WtW = W'*W;
WtX = W'*X;
% Initial error and set of of lambda
err1(1) = max(0,normX2-2*sum(sum(WtX.*H))+sum(sum( WtW.*(H*H'))));
err2(1) = det(WtW); 
lambda = options.lambda * max(1e-6,err1) / (abs( err2 ));
e(1) =  0.5*err1 + 0.5*lambda * err2 ;
 % number of updates of W and H, before the other is updated
if options.display == 1
    disp('Iterations started: '); 
    fprintf('%1.0d ...', 1);
    tic; dispschedtime = 0.1; numdis = 1;  % display parameters
end
% projection model for H
if options.model == 1
    options.proj = 1;
elseif options.model == 2
    options.proj = 2;
elseif options.model == 3
    options.proj = 0;
elseif options.model == 4
    options.proj = 3;    
end
% Main loop 
for i = 2 : options.maxiter
    % *** Update W ***
    for k = 1 : r
        %%%%%%%%%%%%for every wk -start%%%%%%%%%%%%%%
        Ck = null(W(:,[1:k - 1, k + 1:r])');
        gamma = det((W(:,[1:k - 1, k + 1:r])')*W(:,[1:k - 1, k + 1:r]));
        Qk = 0.5 * (H(k, :) * (H(k, :)') * eye(size(X, 1)) + lambda * gamma * Ck * (Ck'));
        Vk = X;
        for t = 1 : r
            if t ~= k
                Vk = Vk - W(:, t) * H(t, :);
            end
        end
        % QQ comment: in the paper, f_k = -Vk * H(k, :)', 
        % which is for (P1) w_k^T Qk w_k + f_k^T x_k. 
        % In the following, we use FGMqpnonnegdet, where 
        % min_{x_i in R^r_+}sum_{i=1}^m ( x_i^T A x_i - 2 c_i^T x_i )is objective function.
        % Thus, here we need to transfer (P1) in the paper 
        % w_k^T Qk w_k + f_k^T x_k 
        % into the form 
        % x_i^T A x_i - 2 c_i^T x_i.
        % Namely, w_k^T Qk w_k -2*(-1/2)f_k^T x_k = w_k^T Qk w_k -2*(-1/2f_k^T) x_k = w_k^T Qk w_k -2*(-1/2*f_k)^T x_k
        % Therefore, we take the fk as -1/2*f_k = 0.5 * Vk * H(k, :)'
        fk = 0.5 * Vk * H(k, :)';
        %QQ add to check whether QK is positive definite
        try 
           chol(Qk);
        catch ME
        disp('Matrix is not symmetric positive definite')
        end
        % Qk of size m\times m, fk of size m\times 1, W(:, k) of size
        % m\times 1
        W(:, k) = FGMqpnonnegforeachcolumnW(Qk,fk,W(:, k),options.inneriter,2);
    %%%%%%%%%%%for every wk--end %%%%%%%%%%%%%%%%%%%%
    end
    % *** Update H ***
    options.init = H; 
    [H,WtW,WtX] = nnls_FPGM(X,W,options);
    err1(i) = max(0, normX2 - 2*sum(sum( WtX.*H ) )  + sum(sum( WtW.*(H*H') ) ) );
    err2(i) = det ( WtW );
    e(i) = 0.5*err1(i) + 0.5*lambda * err2(i);
    if options.display == 1 && (toc >= dispschedtime || i == options.maxiter)
        fprintf('%1.0d ...', i);
        numdis = numdis+1; 
        if mod(numdis,10) == 0
            fprintf('\n');
        end
        tic; 
        dispschedtime = min(20,dispschedtime*1.5); 
    end
    % Tuning lambda to obtain options.target relative error 
    if isfield(options,'target')
        if sqrt(err1(i))/normX > options.target+0.001
            lambda = lambda*0.95;
        elseif sqrt(err1(i))/normX < options.target-0.001
            lambda = lambda*1.05;
        end
    end
end
