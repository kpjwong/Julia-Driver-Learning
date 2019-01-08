function zone2coord(z)
  if z==2
      i = rand(29:45);
      j = rand(107:110);
  elseif z==3
      i = rand(17:28);
      j = rand(60:110);
  elseif z==4
      i = rand(29:36);
      j = rand(39:58);
  elseif z==5
      i = rand(37:45);
      j = rand(43:58);
  elseif z==6
      i = rand(29:36);
      j = rand(15:38);
  elseif z==7
      i = rand(1:17);
      j = rand(97:110);
  elseif z==8
      i = rand(1:8);
      j = rand(1:14);
  elseif z==9
      i = rand(21:28);
      j = rand(15:26);
  elseif z==10
      i = rand(21:28);
      j = rand(35:42);
  elseif z==11
      u = rand();
      i = (u>.5)*rand(9:12)+(u<=.5)*rand(1:12);
      j = (u>.5)*rand(25:28)+(u<=.5)*rand(15:24);
  elseif z==12
      i = rand(9:24);
      j = rand(9:14);
  elseif z==14
      i = rand(1:8);
      j = rand(25:34);
  elseif z==15
      i = rand(1:8);
      j = rand(60:79);
  elseif z==16
      i = rand(29:36);
      j = rand(59:73);
  elseif z==17
      i = rand(37:45);
      j = rand(59:73);
  elseif z==18
      i = rand(29:45);
      j = rand(98:106);
  elseif z==19
      i = rand(13:20);
      j = rand(43:56);
  elseif z==20
      i = rand(9:20);
      j = rand(43:59);
  elseif z==21
      u = rand();
      i = (u>.5)*rand(21:28)+(u<=.5)*rand(13:20);
      j = (u>.5)*rand(55:59)+(u<=.5)*rand(57:59);
  elseif z==22
      u = rand();
      i = (u>.5)*rand(17:20)+(u<.5)*rand(13:20);
      j = (u>.5)*rand(36:42)+(u<.5)*rand(29:35);
  elseif z==23
      u = rand(1:3);
      i = (u==1)*rand(1:8)+(u==2)*rand(9:16)+(u==3)*rand(9:12);
      j = (u==1)*rand(35:39)+(u==2)*rand(36:42)+(u==3)*rand(35:39);
  elseif z==24
      u = rand(1:3);
      i = (u==1)*rand(21:28)+(u==2)*rand(1:8)+(u==3)*rand(9:12);
      j = (u==1)*rand(27:34)+(u==2)*rand(35:39)+(u==3)*rand(29:35);
  elseif z==25
      i = rand(1:8);
      j = rand(50:59);
  elseif z==26
      i = rand(21:28);
      j = rand(43:54);
  elseif z==27
      i = rand(1:8);
      j = rand(40:49);
  elseif z==28
      i = rand(13:20);
      j = rand(15:28);
  elseif z==29
      i = rand(9:17);
      j = rand(78:96);
  elseif z==30
      i = rand(9:17);
      j = rand(60:77);
  elseif z==31
      i = rand(29:45);
      j = rand(87:97);
  elseif z==32
      i = rand(29:45);
      j = rand(74:86);
  elseif z==33
      i = rand(37:45);
      j = rand(15:42);
  elseif z==35
      i = rand(1:8);
      j = rand(80:96);
  else
      i = rand(1:45);
      j = rand(1:110);
  end
  i,j
end
