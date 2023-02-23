function out = Lambdasimilarity(x,y)

out = diag(x' * y ./ (sum(x)'*sum(y)));

end
