function ezcentroid(z)
  ci = Int(5+(z-1-floor((z-1)/5)*5)*9)
  cj = Int(11*ceil(z/5)-5)
  [ci cj]
end
