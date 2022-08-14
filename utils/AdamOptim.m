function [theta_best,loss] = AdamOptim(para,idd_i,optimOptions,f,g)

    T = optimOptions.Niter;
    alpha = optimOptions.alpha;
    beta1 = optimOptions.beta1;%0.9
    beta2 = optimOptions.beta2;%0.999
    epsilon = optimOptions.epsilon;%1e-8
    tol = optimOptions.tol;%1e-6;
    lb = optimOptions.lb;% lower bound
    ub = optimOptions.ub;% upper bound
    loss = zeros(T,1);
    m_tm1 = 0;
    v_tm1 = 0;
    theta_tm1 = para;
    
    theta_best = para;
    loss_best = 1e9;
    loss(1) = norm((f(theta_tm1) - idd_i),'fro');
    for t = 2:T
        % get gradient = jacobian*error
        g_t = 2*g(theta_tm1)'*(f(theta_tm1) - idd_i);
        % Update biased first moment estimate
        m_t = beta1*m_tm1 + (1-beta1)*g_t;
        % Update biased second raw moment estimate
        v_t = beta2*v_tm1 + (1-beta2)*g_t.^2;
        % Compute bias-corrected first moment estimate
        m_t_hat = m_t / (1-beta1^(t-1));
        % Compute bias-corrected second raw moment estimate
        v_t_hat = v_t / (1-beta2^(t-1));
        % Update parameters
        theta_t = theta_tm1 - alpha*m_t_hat./(sqrt(v_t_hat)+epsilon);
        
        % constrain
        theta_t(theta_t < lb) = lb(theta_t < lb);
        theta_t(theta_t > ub) = ub(theta_t > ub);
        
        theta_tm1 = theta_t;
        m_tm1 = m_t;
        v_tm1 = v_t;
            
        idd_pred = f(theta_t);
        loss(t) = norm((idd_pred - idd_i),'fro');

        if mod(t,100) == 0
            alpha = alpha*0.9;
        end
        
        if loss(t) < loss_best
           loss_best = loss(t);
           theta_best = theta_t;
        end
        if (abs(loss(t) - loss(t-1)) < tol)
            loss = loss(1:t);
            break;
        end
    end
    
end