function vec_coord2zone(i,j)
  n = size(i,1)
  z = zeros(n,1)
  for k = 1:n
    if j[k]>=107 && j[k]<=110 && i[k]>=29 && i[k]<=45
        z[k] = 2;
    elseif j[k]>=60 && j[k]<=110 && i[k]>=17 && i[k]<=28
        z[k] = 3;
    elseif j[k]>=39 && j[k]<=58 && i[k]>=29 && i[k]<=36
        z[k] = 4;
    elseif j[k]>=43 && j[k]<=58 && i[k]>=37 && i[k]<=45
        z[k] = 5;
    elseif j[k]>=15 && j[k]<=38 && i[k]>=29 && i[k]<=36
        z[k] = 6;
    elseif j[k]>=97 && j[k]<=110 && i[k]>=1 && i[k]<=17
        z[k] = 7;
    elseif j[k]>=1 && j[k]<=14 && i[k]>=1 && i[k]<=8
        z[k] = 8;
    elseif j[k]>=15 && j[k]<=26 && i[k]>=21 && i[k]<=28
        z[k] = 9;
    elseif j[k]>=35 && j[k]<=42 && i[k]>=21 && i[k]<=28
        z[k] = 10;
    elseif (j[k]>=25 && j[k]<=28 && i[k]>=9 && i[k]<=12) || (j[k]>=15 && j[k]<=24 && i[k]>=1 && i[k]<=12)
        z[k] = 11;
    elseif j[k]>=9 && j[k]<=14 && i[k]>=9 && i[k]<=24
        z[k] = 12;
    elseif j[k]>=25 && j[k]<=34 && i[k]>=1 && i[k]<=8
        z[k] = 14;
    elseif j[k]>= 60 && j[k]<=79 && i[k]>=1 && i[k]<=8
        z[k] = 15;
    elseif j[k]>=59 && j[k]<=73 && i[k]>=29 && i[k]<=36
        z[k] = 16;
    elseif j[k]>=59 && j[k]<=73 && i[k]>=37 && i[k]<=45
        z[k] = 17;
    elseif j[k]>=98 && j[k]<=106 && i[k]>=29 && i[k]<=45
        z[k] = 18;
    elseif j[k]>=43 && j[k]<=56 && i[k]>=13 && i[k]<=20
        z[k] = 19;
    elseif j[k]>=43 && j[k]<=59 && i[k]>=9 && i[k]<=12
        z[k] = 20;
    elseif (j[k]>=55 && j[k]<=59 && i[k]>=21 && i[k]<=28) || (j[k]>=57 && j[k]<=59 && i[k]>=13 && i[k]<=20)
        z[k] = 21;
    elseif (j[k]>=36 && j[k]<=42 && i[k]>=17 && i[k]<=20) || (j[k]>=29 && j[k]<=35 && i[k]>=13 && i[k]<=20)
        z[k] = 22;
    elseif (j[k]>=35 && j[k]<=39 && i[k]>=1 && i[k]<=8) || (j[k]>=36 && j[k]<=42 && i[k]>=9 && i[k]<=16) || (j[k]>=29 && j[k]<=35 && i[k]>=9 && i[k]<=12)
        z[k] = 23;
    elseif (j[k]>=27 && j[k]<=34 && i[k]>=21 && i[k]<=28) || (j[k]>=35 && j[k]<=39 && i[k]>=1 && i[k]<=8) || (j[k]>=29 && j[k]<=35 && i[k]>=9 && j[k]<=12)
        z[k] = 24;
    elseif j[k]>=50 && j[k]<=59 && i[k]>=1 && i[k]<=8
        z[k] = 25;
    elseif j[k]>=43 && j[k]<=54 && i[k]>=21 && i[k]<=28
        z[k] = 26;
    elseif j[k]>=40 && j[k]<=49 && i[k]>=1 && i[k]<=8
        z[k] = 27;
    elseif j[k]>=15 && j[k]<=28 && i[k]>=13 && i[k]<=20
        z[k] = 28;
    elseif j[k]>=78 && j[k]<=96 && i[k]>=9 && i[k]<=17
        z[k] = 29;
    elseif j[k]>=60 && j[k]<=77 && i[k]>=9 && i[k]<=17
        z[k] = 30;
    elseif j[k]>=87 && j[k]<=97 && i[k]>=29 && i[k]<=45
        z[k] = 31;
    elseif j[k]>=74 && j[k]<=86 && i[k]>=29 && i[k]<=45
        z[k] = 32;
    elseif j[k]>=15 && j[k]<=42 && i[k]>=37 && i[k]<=45
        z[k] = 33;
    elseif j[k]>=80 && j[k]<=96 && i[k]>=1 && i[k]<=8
        z[k] = 35;
    else
        z[k] = 0;
    end
  end
  z = convert(Array{Int64,2},z)
end
