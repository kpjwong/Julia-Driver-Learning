function impossible(z)
d = [1]
z = z[1]
  if z<=5
    d = hcat(d,3)
  end
  if z>=46
    d = hcat(d,2)
  end
  if mod(z,5)==1
    d = hcat(d,5)
  end
  if mod(z,5)==0
    d = hcat(d,4)
  end
  if d==[1]
    d = []
  else
    d = d[:,2:end]
  end
  d
end
