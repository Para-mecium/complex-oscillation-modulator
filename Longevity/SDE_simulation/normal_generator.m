function Z = normal_generator(dimension,N, mu, sigma)
    % N: Number of random samples to generate
    % mu, sigma: Mean and standard deviation of the normal distribution

    % Initialize an array to store the generated random numbers
    Z = zeros(dimension, N);

    for i = 1:N
        Z(:,i) = mu+sigma*randn(dimension,1);
    end
end
