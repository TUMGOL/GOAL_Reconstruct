%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
addpath('./utilis/');
close all;
IName = 'lena.png';
%IName = 'august.png';

load('GOAL_OPERATOR');
psz = sqrt(size(Omega,2));

% This is the image we will work on
Ground_truth         = ((imread(['./Images/',IName])));
if size(Ground_truth,3)~=1
    Ground_truth = rgb2gray(Ground_truth);
end
Ground_truth = double(Ground_truth);

range = [0,255];

randn('seed',0);
% Here we set that we will remove 90% of all pixels
missing = 0.9;
% This operator M does the removal for us
M       = Extraction(get_rand(numel(Ground_truth), round((1-missing)*numel(Ground_truth))), size(Ground_truth));
% Remove it by simple matrix vector notation
y       = M*Ground_truth;
% Initialize all missing pixels with something, here we used interpolation
x0      = M'*y;
[Xq,Yq] = meshgrid(1:size(Ground_truth,1),1:size(Ground_truth,2));
X = M*Xq;
Y = M*Yq;
Z = y;
x0 = griddata(X,Y,Z,Xq,Yq);
x0(isnan(x0))=mean(y(:));

l1Smooth     = 1;
lambda0      =  100;
lambda1      =  10;

AT = eye(psz^2)-1/(psz^2);

op = FullOp(Omega*AT, psz, 1, 1, size(Ground_truth));

SR_Instance  =  SparseSolver(y, x0, M, op, l1Smooth, 2, range);

% P-norm used for the recovery process
SR_Instance.p_sm           = 2;

% Maximum number of CG iterations
SR_Instance.max_iterations = 50;

% Lambdas at beginning and end of continuation
SR_Instance.lambda         = lambda0;
SR_Instance.lambda_min     = lambda1;

% Some evaluation things
% These "functions" are called every mod(iteration, SR_Instance.verbose)
% iteration
SR_Instance.side_eval_function{1} =  @(x)psnr(uint8(Ground_truth(:,:,:)),uint8(x(:,:,:)),255);
SR_Instance.verbose               =  1;

% Maximum number of continuation steps
max_out_iter  = 3;

% Continuation parameter for the lagrange multiplier lambda
mul_la        = (SR_Instance.lambda_min/SR_Instance.lambda)^(1/(max_out_iter-1));

% The conjugate gradient update formula to be applied
CG_method     = 'HH';

SR_Instance.maxLineSearchSteps = 50;

for cr = 1:max_out_iter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Here is where the magic takes place
    SR_Instance.t = max(SR_Instance.t,10);
    conjugate_gradient(SR_Instance, CG_method, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SR_Instance.lambda   = mul_la * SR_Instance.lambda;
    fprintf('PSNR PSNR PSNR\n')
    fprintf('PSNR PSNR PSNR\n')
    psnr((Ground_truth(1:end-0,1:end-0)),double(uint8((SR_Instance.x(1:end,1:end)))),255);
    fprintf('PSNR PSNR PSNR\n')
    fprintf('PSNR PSNR PSNR\n')
    
end

Rec = SR_Instance.x;
Rec(1:2,:) = repmat(SR_Instance.x(3,:),2,1);
Rec(end-1:end,:) = repmat(SR_Instance.x(end-2,:),2,1);
Rec(:,1:2) = repmat(SR_Instance.x(:,3),1,2);
Rec(:,end-1:end) = repmat(SR_Instance.x(:,end-2),1,2);

figure(100)
% Original Image
subplot(3,1,1)
imshow(Ground_truth/max(abs(Ground_truth(:))))
xlabel('Orignial Image')

% Starting Point
subplot(3,1,2)
imshow(x0/max(abs(x0(:))))
xlabel('Initial Guess Image X0')

% Reconstruction
subplot(3,1,3)
imshow(Rec/max(abs(Rec(:))))
xlabel('Reconstruction Result')







%


