function eztozone(zz,dir)
  zz = zz[1]
  dir = dir[1]
  if dir==1
    z = zz
  elseif dir==2
    zz>=46 ? z = NaN : z = zz+5
  elseif dir==3
    zz<=5 ? z = NaN : z = zz-5
  elseif dir==4
    mod(zz,5)==0 ? z = NaN : z = zz+1
  elseif dir==5
    mod(zz,5)==1 ? z = NaN : z = zz-1
  end
  if isnan(z)==0
    z = Int(z)
  end
  z
end
