function ezvec_coord2zone(a)
  i = a[:,1]
  j = a[:,2]
  n = size(a,1)
  z = zeros(n,1)
  for k = 1:n
    z[k] = 5*(ceil(j[k]/11)-1)+ceil(i[k]/9)
  end
  z = convert(Array{Int64,2},z)
  z
end
