function B = shr(A, k)

if k < 0
  B = shl(A, -k);
  return
end

B = A(:, 1 : size(A, 2) - k);
B = cat(2, zeros(size(A, 1), size(A, 2) - size(B, 2)), B);

end
