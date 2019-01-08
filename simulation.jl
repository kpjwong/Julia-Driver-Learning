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
r = reshape([],0,1)
W = []
S = []
t = 1
h = 0
c = .0057
u = 15
zm = reshape([],0,1)
cx = rand(1:45,Nc[h+1]-1,1)
cy = rand(1:110,Nc[h+1]-1,1)
C = hcat(cx,cy,vec_coord2zone(cx,cy),zeros(Nc[h+1]-1,2))
CC = reshape([],0,7)
CS = reshape([],0,7)
D = reshape([],0,4)
qx = rand(1:45)
qy = rand(1:110)
P = reshape(readdlm("C:\\Users\\jerem\\Documents\\MATLAB\\P.txt", ','),45,110,24)
dec1 = map(i->rand()<P[cx[i],cy[i],h+1],1:size(C,1))
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
DD = vcat(map((i,j)->repmat([i j vec_coord2zone(i,j) 0],l[j,i],1),repeat(1:45,inner=110),repeat(1:110,outer=45))...)
D = vcat(D,DD)
Q_payoff = 0
old_Q_payoff = 0
tau = 0
qS = reshape([],0,3)
qs = 0
qa = 0
PI = []

iteration = 1
match = reshape([],0,2)
agg_zone_match = reshape([],0,3)
agg_r = reshape([],0,2)
pct = []
for iteration = 1:300
  for h = 0:23
    m = 0
    r = []
    if Nc[h+1] > size(C,1)+size(CC,1)
      Nc_diff = Nc[h+1]-size(C,1)-size(CC,1)
      cx = rand(1:45,Nc_diff-1,1)
      cy = rand(1:110,Nc_diff-1,1)
      new_dec1 = map(i->rand()<P[cx[i],cy[i],h+1],1:Nc_diff-1)
      new_dec2 = (-1).^rand(1:2,Nc_diff-1,1)
      new_dec2[intersect(find(cx.==1),find(new_dec1.==1))] = 1
      new_dec2[intersect(find(cx.==45),find(new_dec1.==1))] = -1
      new_dec2[intersect(find(cy.==1),find(new_dec1.==0))] = 1
      new_dec2[intersect(find(cy.==110),find(new_dec1.==0))] = -1
      C_diff = hcat(cx,cy,vec_coord2zone(cx,cy),zeros(Nc_diff-1,2),new_dec1,new_dec2)
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
        PP = cumsum(hcat(PP11[Int(cx[n]),Int(cy[n]),Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP12[Int(cx[n]),Int(cy[n]),Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP21[Int(cx[n]),Int(cy[n]),Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1],PP22[Int(cx[n]),Int(cy[n]),Int(2-old_dec1[n]),Int(1.5-.5*old_dec2[n]),h+1]),2)
        draw = sum(rn.>PP)+1
        dec[n,:] = d_ref[draw,:]
      end
      dec1 = dec[:,1]
      dec2 = dec[:,2]
      cx = convert(Array{Int64,1},cx+dec1.*dec2)
      cy = convert(Array{Int64,1},cy+(1.-dec1).*dec2)
      C = hcat(cx,cy,vec_coord2zone(cx,cy),s,payoff,dec1,dec2)
      qz = ezvec_coord2zone([qx qy])
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
          temp_V = Q_stay[qz,h+1]
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
            while sum(impossible(qz).==qa)~=0
              qa = rand(1:5)
            end
          else
            Q = vcat(Q_stay[qz,h+1],zeros(4,1))
            for dir = 2:5
              if sum(impossible(qz).==dir)==0
                centroid = ezcentroid(eztozone(qz,dir))
                ci = centroid[1]
                cj = centroid[2]
                dd = abs(qx-ci)+abs(qy-cj)
                Q[dir] = -c*dd+beta^dd*V[eztozone(qz,dir)]
              else
                Q[dir] = NaN
              end
            end
            mxQ, qa = findmax(Q)
          end
          if qa==1
            cruise_route = cruise(qx,qy,cruise_t,v[h+1])
            tau = -size(cruise_route,1)
          else
            centroid = ezcentroid(eztozone(qz[1],qa[1]))
            ci = centroid[1]
            cj = centroid[2]
            cruise_route = ezroute(qx,qy,ci,cj,qa)
            tau = -size(cruise_route,1)
          end
        end
        if sum(Int[[qx qy] == D[i,1:2]' for i=1:size(D,1)])!=0
          comp_c = sum(Int[[qx qy] == C[i,1:2]' for i=1:size(C,1)],1)
          comp_d = sum(Int[[qx qy] == D[i,1:2]' for i=1:size(D,1)],1)
          if comp_c+1 > comp_d
            comp_draw = rand(comp_c,1)
            qdraw = sum(comp_draw.<rand(),1)+1
          end
          if comp_c+1 <= comp_d || qdraw > comp_c-comp_d
            ind_d = findin(Int[[qx qy] == D[i,1:2]' for i=1:size(D,1)],1)
            D = D[deleteat!(collect(1:size(D,1)),ind_d),:]
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
              temp_V = Q_stay[qz,h+1];
              for dir = 2:5
                if sum(impossible(qz).==dir)==0
                  temp_dist = 13-(dir-mod(dir,2));
                  temp_V = vcat(temp_V,-c*temp_dist+beta^temp_dist*V[eztozone(qz,dir),h+1]);
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
        ~,ind_c,ind_d = intersectML(C[:,1:2],D[:,1:2])
        while !isempty(ind_c) && !isempty(ind_d)
          W = vcat(W,D[ind_d,3],D[ind_d,4]/T[h+1]*60,h*ones(size(ind_d,1),1))
          S = vcat(S,C[ind_c,3],C[ind_c,4]/T[h+1]*60,h*ones(size(ind_c,1),1))
          type_draw = rand(size(ind_c))
          trip_type = map(x->sum(x.>trip_type_draw[h+1,:]),type_draw)
          temp_CC = [C[ind_c,1:3] C[ind_c,5]]
          cc_xy = zeros(size(temp_CC,1),2)
          cc_z_tlc = zeros(size(trip_type,1),1)
          ng_zone_draw = min(map(x->sum(x.>tp_ng'),rand(sum(trip_type.==1),1)).+1,33)
          nm_zone_draw = min(map(x->sum(x.>tp_nm'),rand(sum(trip_type.==2),1)).+1,33)
          cc_z_tlc[trip_type.==1] = zone_id[ng_zone_draw]
          cc_z_tlc[trip_type.==2] = zone_id[nm_zone_draw]
          temp_cc_z3 = temp_CC[trip_type.==3,3]
          cc_z_tlc3 = map(n->zone_id[min(sum(rand().>tp[:,find(i->zone_id[i]==temp_cc_z3[n],1:33),h+1])+1,33)],1:sum(trip_type.==3))
          cc_z_tlc[trip_type.==3] = cc_z_tlc3
          cc_xy = map(x->zone2coord(x),cc_z_tlc)
          cc_xy = vcat(map(i->[cc_xy[i][1] cc_xy[i][2]],1:size(cc_xy,1))...)
          timer = zeros(size(temp_CC,1),1)
          timer[trip_type.==1] = max(20*ceil(min(rand(LogNormal(mu_nongrid[h+1],sig_nongrid[h+1],sum(trip_type.==1),1)),20)),4)
          timer[trip_type.==2] = max(20*ceil(min(rand(LogNormal(mu_nonman[h+1],sig_nonman[h+1],sum(trip_type.==2),1)),50)),4)
          timer[trip_type.==3] = max(abs(cc_xy[trip_type.==3,1]-temp_CC[trip_type.==3,1])+abs(cc_xy[trip_type.==3,2]-temp_CC[trip_type.==3,2]),4)
          payoff = temp_CC[:,4]+2.5+.125*timer
          cc_dec1 = map(i->Int(rand()<P[cc_xy[i,1],cc_xy[i,2],h+1]),1:size(cc_xy,1))
          cc_dec2 = (-1).^rand(1:2,size(cc_xy,1),1)
          cc_dec2[intersect(find(cc_xy[:,1].==1),find(cc_dec1.==1))] = 1
          cc_dec2[intersect(find(cc_xy[:,1]==45),find(cc_dec1.==1))] = -1
          cc_dec2[intersect(find(cc_xy[:,2]==1),find(cc_dec1.==0))] = 1
          cc_dec2[intersect(find(cc_xy[:,2]==110),find(cc_dec1.==0))] = -1
          CC = vcat(CC,[cc_xy cc_z_tlc timer payoff cc_dec1 cc_dec2])
          zm = vcat(zm,C[ind_c,3])
          C = C[deleteat!(collect(1:size(C,1)),sort(ind_c)),:]
          D = D[deleteat!(collect(1:size(D,1)),sort(ind_d)),:]
          m = m+size(ind_c,1)
          ~,ind_c,ind_d = intersectML(C,D)
        end
        if !isempty(CC)
          CC = CC-repmat([0 0 0 1 c 0 0],size(CC,1),1)
          timer = CC[:,4]
          CS = CC[timer.==0,:]
          CC = CC[deleteat!(collect(1:size(CC,1)),sort(find(timer.==0))),:]
        end
        C = vcat(C,CS)
        CS = reshape([],0,7)
        r = vcat(r,size(CC,1)/Nc[h+1])
        l = poissrnd(1.3*Grid(lambda[:,h+1])/T[h+1])
        DD = vcat(map((i,j)->repmat([i j vec_coord2zone(i,j) 0],l[j,i],1),repeat(1:45,inner=110),repeat(1:110,outer=45))...)
        D = vcat(D,DD)
        t = t+1
      end
      agg_r = vcat(agg_r,[mean(r),h])
      match = vcat(match,[m h])
      agg_zone_match = vcat(agg_zone_match,map(n->[zone_id[n] sum(zm.==zone_id[n]) h],1:33))
      C[:,5] = zeros(size(C,1),1)
      CC[:,5] = zeros(size(CC,1),1)
      h = h+1
    end
  end
