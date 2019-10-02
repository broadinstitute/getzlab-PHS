function q = fdr_jh(p, n_extra)

if nargin == 1, n_extra = 0; end 
if min(size(p)) ~= 1, error('Only vectors supported at this time.'); end

[sp, ord] = sort(p);
q = [sp*(length(p) + n_extra)./(1:length(p))'; 1];

q(q > 1) = 1;

%ensure values are monotonically decreasing
for i = length(p):-1:1,
  q(i) = min(q([i i + 1]));
end
q(end) = [];

q(ord) = q;
