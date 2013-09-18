function conjugate_gradient(Sparse_Problem, CG_method, Disp_mat, Ground_truth)

Sparse_Problem.init();

%aviobj = avifile('CS_rand.avi','compression','none','fps',10);
if nargin == 4
    figure(100)
    clf
    set(gcf,'Position',[50 50 620 780])
    %subplot(2,2,[3,4])
    
    subplot(2,1,2)
    hold on
    grid on
    XX=Disp_mat*Sparse_Problem.x;
    axis([1,Sparse_Problem.max_iterations,0,mean(mean(abs(XX-Ground_truth)))]);
    ylabel(['Mean Absolute Disparity Error'],'FontWeight','bold','FontSize',12);
    xlabel('Iteration','FontWeight','bold','FontSize',12);
end

%screen_size = get(0, 'ScreenSize');
%f1 = figure(1234);
%set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );



while(1)
    % First we perform the linesearch, where we find the step size t and
    % update the current solution to x
    [Sparse_Problem,worked] = Sparse_Problem.linesearch();
    if ~worked
        return;
    end

    % Compute the gradient at the next point which we have found using the
    % linesearch above
    g1 = Sparse_Problem.gradient();
    g0 = Sparse_Problem.g;
    % Now we update the search direction using any conjugate gradient
    % method
    cg_beta = 1;
    switch CG_method
        case 'HS'
            yk       =  g1 - g0;
            cg_beta  = (real(g1(:)'*yk(:))/real(Sparse_Problem.dx(:)'*yk(:)+eps));
        case 'PR'
            cg_beta  = max(0,(g1(:)'*(g1(:)-g0(:)))/(g0(:)'*g0(:)));
        case 'MK'
            ghs      =    g1 - g0;
            cg_beta  =  -(g1(:)'*ghs(:))/(Sparse_Problem.dx(:)'*g0(:));
        case 'DY'
            yk      =   g1-g0;
            cg_beta  =   (g1(:)'*g1(:))/(Sparse_Problem.dx(:)'*yk(:));
        case 'FR'
            cg_beta  = (g1(:)'*g1(:))/(g0(:)'*g0(:));
        case 'HZ'
            cg_beta  = -(g1(:)'*g1(:))/(Sparse_Problem.dx(:)'*g0(:));
        case 'HH'
            yk       =  g1 - g0;
            denom = real(Sparse_Problem.dx(:)'*yk(:)+eps);
            cg_beta     = (g1(:)'*yk(:))/denom;
            cg_beta_dy  = (g1(:)'*g1(:))/denom;
            cg_beta     = max(0,min(cg_beta_dy, cg_beta));
        case 'SD'
            cg_beta = 0;
    end

    % Update this for the next iteration
    [Sparse_Problem, stop] = Sparse_Problem.next_iteration(cg_beta, g1);
    if stop
        return;
    end
    
    if Disp_mat == 1
        figure(1234)
        if size(Sparse_Problem.x,3) > 3
            slices = 1;%size(Sparse_Problem.x,3);
            for i=1:slices
                subplot(2,slices,i)
                normi = max(max(Sparse_Problem.x(:,:,i)));
                imshow(Sparse_Problem.x(:,:,i)/normi);
                drawnow
                subplot(2,slices,i+slices)
                normi = max(max(Sparse_Problem.y(:,:,i)));
                imshow(Sparse_Problem.y(:,:,i)/normi);
                drawnow
            end
            %imshow((Sparse_Problem.x(:,:,1))./max(abs(vec(Sparse_Problem.x(:,:,1)))));
        else
            subplot(1,2,1)
            imshow(uint8(Sparse_Problem.x(:,:,:)));
            title('Reconstruction Result s')
            subplot(1,2,2)
            imshow(uint8(Sparse_Problem.y));
            title('Measurement y')
            %imshow((Sparse_Problem.x(:,:,1))./max(abs(vec(Sparse_Problem.x(:,:,1)))));
        end
    end
    if nargin == 4
        capture_video;
    end
    
end

end