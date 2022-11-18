startlam = 6;
n_iter   = 20;
tolerance = 1e-1;
lam = startlam;
grad_lam    = 2*lam;
fprintf('gradlam: %f\n',grad_lam);
learn_rate  = 0.4;
%objective func: v^2

[x,y] = grad_descent_lam(grad_lam,...
    startlam,learn_rate,n_iter,tolerance)

function [nbiterations, converged_lam] = grad_descent_lam(grad_lam,...
    lam,learn_rate,n_iter,tolerance)
    nbiterationslam = 1;
    
    for i = 1:n_iter
        grad_lam = 2*lam;%change here
        diff = -learn_rate*grad_lam;
        
        if (abs(diff)<= tolerance)
            converged_lam = lam;
            nbiterationslam = nbiterationslam+1;
            nbiterations =  nbiterationslam;
            disp('yea')
            break;
        else
            converged_lam = 0;
            nbiterationslam = nbiterationslam+1;
            nbiterations =  nbiterationslam;
            lam = lam + diff;
            fprintf('lam: %f\n',lam);
        end
    end
    
end