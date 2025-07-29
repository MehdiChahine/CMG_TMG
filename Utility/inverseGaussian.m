function x = inverseGaussian(mu, lambda)
    % Generate a standard normal random variable
    v = randn;
    y = v^2;
    
    % Compute the candidate value x
    x = mu + (mu^2 * y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu^2 * y^2);
    

    % Accept or transform x
    if rand <= mu / (mu + x)
        return;
    else
        x = (mu^2) / x;
    end
end
