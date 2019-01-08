using QL
using Distributions
param = readdlm("parameters.txt", ',')
v = param[:,1]
p_nongrid = param[:,4]
p_nonman = param[:,5]
trip_type_draw = hcat(zeros(24,1),p_nongrid,p_nongrid+p_nonman)
mu_nongrid = param[:,6]
sig_nongrid = param[:,7]
mu_nonman = param[:,8]
sig_nonman = param[:,9]
lambda = readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\lambda.txt", ',')
zone_id = [2;3;4;5;6;7;8;9;10;11;12;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;35;0]
Nc = convert(Array{Int64,2},floor(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\Nc.txt", ',')))
tp = cumsum(reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\tpp.txt", ','),33,33,24),1)
tp2 = readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\tp_ngnm.txt", ',')
tp_ng = cumsum(tp2[:,1]/100)
tp_nm = cumsum(tp2[:,2]/100)
alpha = .1
beta = .99
cruise_t = 5
epsilon = .01
PP11 = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\PP11.txt", ','),45,110,2,2,24)
PP12 = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\PP12.txt", ','),45,110,2,2,24)
PP21 = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\PP21.txt", ','),45,110,2,2,24)
PP22 = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\PP22.txt", ','),45,110,2,2,24)
d_ref = [1 1;1 -1;0 1;0 -1]
Q_stay = readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\Q_stay.txt", ',')
V = zeros(50,24)
a = zeros(45,110,24)
direction = direction = ['0' 'N' 'S' 'E' 'W']
T = ceil(20*v)
r = []
W = []
S = []
t = 1
h = 0
c = .0057
u = 15

cx = rand(1:45,Nc[h+1]-1,1)
cy = rand(1:110,Nc[h+1]-1,1)
C = hcat(cx,cy,vec_coord2zone([cx cy]),zeros(Nc[h+1]-1,2))
CC = []
CS = []
D = [1 1 1]
qx = rand(1:45)
qy = rand(1:110)
P = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\P.txt", ','),45,110,24)
dec1 = zeros(size(C,1),1)
for i = 1:size(C,1)
  dec1[i] = rand() < P[cx[i],cy[i],h+1]
end
dec2 = (-1).^rand(1:2,size(C,1),1)
dec2[intersect(find(cx.==1),find(dec1.==1))] = 1
dec2[intersect(find(cx.==45),find(dec1.==1))] = -1
dec2[intersect(find(cy.==1),find(dec1.==0))] = 1
dec2[intersect(find(cy.==110),find(dec1.==0))] = -1
C = hcat(C,dec1,dec2)
l = zeros(110,45)
poissrnd(l) = rand(Poisson(l))
@vectorize_1arg Number poissrnd
while l==zeros(110,45)
  l = poissrnd(1.3*Grid(lambda[:,h+1])/T[h+1])
end
for j = 1:110
  for i = 1:45
    D = vcat(D,repmat(hcat(i,j,vec_coord2zone(hcat(i,j)),0),l[j,i],1))
  end
end
D = D[2:end,:]
Q_payoff = 0
old_Q_payoff = 0
tau = 0
qS = reshape([],0,3)
qs = 0
qa = 0
PI = []

iteration = 1
match = []
zone_match = hcat(zone_id,zeros(33,2))
agg_zone_match = []
agg_r = []
pct = []

    m = 0
    r = []
    if Nc[h+1] > size(C,1)+size(CC,1)
      Nc_diff = Nc[h+1]-size(C,1)-size(CC,1)
      cx = rand(1:45,Nc_diff-1,1)
      cy = rand(1:110,Nc_diff-1,1)
      new_dec1 = zeros(Nc_diff-1,1)
      for i = 1:Nc_diff-1
        new_dec1[i] = rand() < P[cx[i],cy[i],h+1]
      end
      new_dec2 = (-1).^rand(1:2,Nc_diff-1,1)
      new_dec2[intersect(find(cx.==1),find(dec1.==1))] = 1
      new_dec2[intersect(find(cx.==45),find(dec1.==1))] = -1
      new_dec2[intersect(find(cy.==1),find(dec1.==0))] = 1
      new_dec2[intersect(find(cy.==110),find(dec1.==0))] = -1
      C_diff = hcat(cx,cy,vec_coord2zone(hcat(cx,cy)),zeros(Nc_diff-1,2),new_dec1,new_dec2)
      C = vcat(C,C_diff)
    elseif Nc[h+1] < size(C,1)+size(CC,1)
      Nc_diff = size(C,1)+size(CC,1)-Nc[h+1]
      C = C[Nc_diff:end,:]
    end
    t = 1
    while t < T[h+1]
      D = D+repmat([0 0 0 1],size(D,1),1)
      C = C+repmat([0 0 0 1 -c 0 0],size(C,1),1)
      cx = C[:,1]
      cy = C[:,2]
      s = C[:,4]
      payoff = C[:,5]
      old_dec1 = C[:,6]
      old_dec2 = C[:,7]
      dec = zeros(size(C,1),2)
      for n = 1:size(C,1)
        rn = rand()
        while rn >= .99999
          rn = rand()
        end
        PP = cumsum(hcat(PP11[cx[n],cy[n],Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP12[cx[n],cy[n],Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP21[cx[n],cy[n],Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP22[cx[n],cy[n],Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1]),2)
        draw = sum(rn.>PP)+1
        dec[n,:] = d_ref[draw,:]
      end
      dec1 = dec[:,1]
      dec2 = dec[:,2]
      cx = cx+dec1.*dec2
      cy = cy+(1.-dec1).*dec2
      C = hcat(cx,cy,vec_coord2zone(hcat(cx,cy)),s,payoff,dec1,dec2)
      qz = ezvec_coord2zone(hcat(qx,qy))
      if tau<=0
        if tau<0
          qx = cruise_route[end+tau+1,1]
          qy = cruise_route[end+tau+1,2]
          tau = tau+1
          qs = qs+1
          Q_payoff = Q_payoff-c
        elseif tau==0
          if qa==1
            Q_stay[qz,h+1] = alpha*(Q_payoff-old_Q_payoff)+(1-alpha)*Q_stay[qx,h+1]
          end
          temp_V = [Q_stay[qz,h+1]]
          for dir = 2:5
            if sum(impossible(qz).==dir)==0
              temp_dist = 13-(dir-mod(dir,2));
              temp_V = vcat(temp_V,-c*temp_dist+beta^temp_dist*V[eztozone(qz,dir),h+1])
            end
          end
          V[qz,h+1] = maximum(temp_V)
          qs = qs+1
          old_Q_payoff = Q_payoff
          Q_payoff = Q_payoff-c
          if rand() < epsilon
            qa = rand(1:5)
            while sum(impossible(qz)==qa)~=0
              qa = rand(1:5)
            end
          else
            Q = vcat(Q_stay[qz,h+1],zeros(4,1))
            for dir = 2:5
              if sum(impossible(qz)==dir)==0
                centroid = ezcentroid(eztozone(qz,dir))
                ci = centroid[1]
                cj = centroid[2]
                dd = abs(qx-ci)+abs(qy-cj)
                Q[dir] = -c*dd+beta^dd*V[eztozone(qz,dir)]
              end
            end
            mxQ, qa = findmax(Q,1)
          end
          if qa==1
            cruise_route = cruise(qx,qy,cruise_t,v[h+1])
            tau = -size(cruise_route,1)
          else
            centroid = ezcentroid(eztozone(qz,qa))
            ci = centroid[1]
            cj = centroid[2]
            cruise_route = ezroute(qx,qy,ci,cj,qa)
            tau = -size(cruise_route,1)
          end
        end
        if sum(Int[[qx qy] == D[i,:] for i=1:size(D,1)])!=0
          comp_c = sum(Int[[qx qy] == C[i,1:2] for i=1:size(C,1)],1)
          comp_d = sum(Int[[qx qy] == D[i,1:2] for i=1:size(D,1)],1)
          if comp_c+1 > comp_d
            comp_draw = rand(comp_c,1)
            qdraw = sum(comp_draw.<rand(),1)+1
          end
          if comp_c+1 <= comp_d || qdraw > comp_c-comp_d
            ind_d = findin(Int[[qx qy] == D[i,1:2] for i=1:size(D,1)],1)
            D = D[deletat!(collect(1:size(D,1)),ind_d),:]
            qS = vcat(qS,[qs*3/v[h+1] ezvec_coord2zone([qx qy]) h])
            trip_type = sum(rand().>trip_type_draw[h+1,:])
              if trip_type==1
                tau = 2*ceil(20*min(rand(LogNormal(mu_nongrid[h+1],sig_nongrid[h+1])),20))
                tlc_qz = zone_id[sum(rand().>tp_ng)+1]
                q_xy = zone2coord(tlc_qz)
                qx = q_xy[1]
                qy = q_xy[2]
              elseif trip_type==2
                tau = 2*ceil(20*min(rand(LogNormal(mu_nonman[h+1],sig_nonman[h+1])),20))
                tlc_qz = zone_id[sum(rand().>tp_nm)+1]
                q_xy = zone2coord(tlc_qz)
                qx = q_xy[1]
                qy = q_xy[2]
              elseif trip_type==3
                tlc_qz = zone_id[sum(rand().>cumsum(tp[:,findin(Int[cood2zone([qx qy]) == zone_id[i,:] for i=1:size(zone_id,1)],1)]))+1]
                old_qx = qx
                old_qy = qy
                q_xy = zone2coord(tlc_qz)
                qx = q_xy[1]
                qy = q_xy[2]
                tau = abs(qx-old_qx)+abs(qy-old_qy)
              end
              Q_payoff = Q_payoff+2.5+(.125-c)*tau+beta^tau*V[ezvec_coord2zone([qx qy]),h+1]
              Q_stay[qz,h+1] = alpha*(Q_payoff-old_Q_payoff)+(1-alpha)*Q_stay[qz,h+1]
              temp_V = [Q_stay(qz,h+1)];
              for dir = 2:5
                if sum(impossible(qz).==dir)==0
                  temp_dist = 13-(dir-mod(dir,2));
                  temp_V = vcat(temp_V,-c*temp_dist+beta^temp_dist*V(eztozone(qz,dir),h+1));
                end
              end
              V[qz,h+1] = maximum(temp_V,1);
            end
          end
        end
        if tau > 0
          tau = tau-1
          qs = 0
        end
