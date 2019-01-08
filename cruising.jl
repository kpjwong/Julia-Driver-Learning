function ezvec_coord2zone(a)
  i = a[:,1]
  j = a[:,2]
  n = size(a,1)
  z = zeros(n,1)
  for k = 1:n
    z[k] = 5*(ceil(j[k]/11)-1)+ceil(i[k]/9)
  end
  convert(Array{Int64,2},z)
end

function ezroute(i1,j1,i2,j2,dir)
  r_start = reshape([],0,2)
  r = reshape([],0,2)
  r_end = reshape([],0,2)
  or_i1 = i1
  or_i2 = i2
  if abs(i1-i2) < 4 && j1==j2
    r = reshape([],0,2)
    theta = sign(i1-i2)
    k = abs(i1-i2)
    for n = 1:k
      r = vcat(r,[i1+n*theta j1])
    end
  else
    if mod(i1-1,4)!=0
      if i1-i2>0
        K1 = mod(i1-1,4)
        theta1 = -1
      else
        K1 = 4-mod(i1-1,4)
        theta1 = 1
      end
      for k = 1:K1
        r_start = vcat(r_start,[i1+k*theta1 j1])
      end
    else
      K1 = 0
      theta1 = 0
    end
    i1 = or_i1+K1*theta1
    if mod(i2-1,4)!=0
      if i1>i2
        K2 = 4-mod(i2-1,4)
        theta2 = 1
      else
        K2 = mod(i2-1,4)
        theta2 = -1
      end
      for k = K2-1:-1:0
        r_end = vcat(r_end,[i2+k*theta2 j2])
      end
    else
      K2 = 0
      theta2 = 0
    end
    i2 = or_i2+K2*theta2
    nx = Int((i2-i1)/4)
    theta = sign(nx)
    nx = abs(nx)
    ny = Int(j2-j1)
    phi = sign(ny)
    ny = abs(ny)
    rx = repmat([1 0],4,1)
    ry = [0 1]
    if dir==2 || dir==3
      ref = vcat(repmat([2],ny,1),ones(Int,nx,1))
    elseif dir==4 || dir==5
      ref = vcat(ones(Int,nx,1),repmat([2],ny,1))
    end
    draw = ref
    for m = 1:size(draw,1)
      if draw[m]==1
        r = vcat(r,theta*rx)
      elseif draw[m]==2
        r = vcat(r,phi*ry)
      end
    end

    if isempty(r)==0
      r = cumsum(r,1)
      r = r+repmat([i1 j1],size(r,1),1)
    end
    r = vcat(r_start,r,r_end)
    if isempty(r)==0 & r[1,:]'==[or_i1 j1]
      r = r[2:end,:]
    end
  end
  r
end

function cruise(i,j,t,v)
  length = ceil(t*v/3)
  z = ezvec_coord2zone(hcat(i,j))
  r = reshape([],0,2)
  l = 0
  dir = rand(2:5)
  D = [0 1;0 -1;1 0;-1 0]
  while l<length
    indx = Int(z[1]-floor((z[1]-1)/5)*5)
    indy = Int(ceil(z[1]/5))
    i1 = rand((indx-1)*9+1:indx*9)
    j1 = rand((indy-1)*11+1:indy*11)
    route = ezroute(i,j,i1,j1,dir[1])
    counter = 0
    if l>1
      dd = r[end,:]-r[end-1,:]
      dir = findin(Int[dd == D[i,:] for i=1:size(D,1)],1)+1
      while ezvec_coord2zone(hcat(i,j)+D[dir[1]-1,:]')!=z
        dir = rand(2:5)
      end
    else
      dir = rand(2:5)
    end
    if l==length-1
      route = r[end,:]'+dd'
      while ezvec_coord2zone(route)!=repmat(z,size(route,1),1)
        route = route[end,:]'+D[rand(1:4),:]'
      end
    else
      if size(route,1)>length-l
        while size(route,1)!=length-l
          i1 = rand((indx-1)*9+1:indx*9)
          j1 = rand((indy-1)*11+1:indy*11)
          route = ezroute(i,j,i1,j1,dir[1])
          counter = counter+1
        end
      end
    end
    r = vcat(r,route)
    i = i1
    j = j1
    l = size(r,1)
  end
  r
end
