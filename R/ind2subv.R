ind2subv <- function(siz, ndx) stop("Needs translation")
# function sub = ind2subv(siz, ndx)
# % IND2SUBV Like the built-in ind2sub, but returns the answer as a row vector.
# % sub = ind2subv(siz, ndx)
# %
# % siz and ndx can be row or column vectors.
# % sub will be of size length(ndx) * length(siz).
# %
# % Example
# % ind2subv([2 2 2], 1:8) returns
# %  [1 1 1
# %   2 1 1
# %   ...
# %   2 2 2]
# % That is, the leftmost digit toggle fastest.
# %
# % See also SUBV2IND

# n = length(siz);

# if n==0
#   sub = ndx;
#   return;
# end

# if all(siz==2)
#   sub = dec2bitv(ndx-1, n);
#   sub = sub(:,n:-1:1)+1;
#   return;
# end

# cp = [1 cumprod(siz(:)')];
# ndx = ndx(:) - 1;
# sub = zeros(length(ndx), n);
# for i = n:-1:1 % i'th digit
#   sub(:,i) = floor(ndx/cp(i))+1;
#   ndx = rem(ndx,cp(i));
# end
