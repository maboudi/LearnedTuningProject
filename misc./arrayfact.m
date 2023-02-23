
function y = arrayfact(x)

% This function calculate elementwise factorial of array x

x_1D = x(:);
for i = 1 : length(x_1D)
    y_1D(i) = factorial(x_1D(i));
end

y = reshape(y_1D, size(x));
end


