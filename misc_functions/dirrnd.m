function d = dirrnd(a, n)

if ~exist('n', 'var'), n = 1; end
if all(size(a) ~= 1), error('Only single parameters supported at this time.'); end

gr = NaN(n, max(size(a)));
for i = 1:n,
  gr(i, :) = gamrnd(a, 1);
end

d = bsxfun(@rdivide, gr, sum(gr, 2));
