from MarkovModels import GenerateQ_11rates
from numpy import array
from scipy.optimize import minimize
from MarkovModels import GetStationaryDistribution
from MarkovModels import GenerateUnrestRateMatrix
from random import uniform
# def GenerateUnrestRateMatrix(pi):

pi = [0.25, 0.4, 0.2, 0.15]
Q = GenerateUnrestRateMatrix(pi)
pi_stat = GetStationaryDistribution(Q)
print(pi_stat)

# def evolve_pi(pi_current,change_gc,sign):
change_gc = 0.1
sign = -1
delta_gc = change_gc * sign
delta_at = change_gc * sign * -1
if sign == -1:
    delta_g = uniform(delta_gc,0)
    delta_a = uniform(0,delta_at)    
else:
    assert(sign == 1)
    delta_g = uniform(0,delta_gc)    
    delta_a = uniform(delta_at,0)
delta_c = delta_gc - delta_g
delta_t = delta_at - delta_a
pi_new = pi
pi_new[0] += delta_a
pi_new[1] += delta_t
pi_new[2] += delta_g
pi_new[3] += delta_c    

    

# pi = [0.25, 0.4, 0.2, 0.15]
# def penalty(rates_11):
#     A2C, A2G, A2T, C2A, C2G, C2T, G2A, G2C, G2T, T2A, T2C = rates_11
#     T2G = 1
#     Q = array([[-(A2C+A2G+A2T),A2C,A2G,A2T],
#                 [C2A,-(C2A+C2G+C2T),C2G,C2T],
#                 [G2A,G2C,-(G2A+G2C+G2T),G2T],
#                 [T2A,T2C,T2G,-(T2A+T2C+T2G)]])
#     pi_stat = GetStationaryDistribution(Q)
#     penalty = 0
#     for i in range(4):        
#         penalty += pow(pi[i] -pi_stat[i],2)
#     return (penalty)
        
# rates_11 = [1]*11
# opt_res = minimize(penalty,rates_11,method='Nelder-Mead')
# # print ("Penalty is ", penalty(opt_res.x))
# # print(opt_res.x)
# for i in range(11):
#     rates_11[i] = opt_res.x[i]
# # rates_11.append(1)
# Q = GenerateQ_11rates(rates_11)
# # pi_stat = GetStationaryDistribution(Q)
# # print(pi_stat[0], pi_stat[1], pi_stat[2], pi_stat[3])
# rootProb = GetStationaryDistribution(Q)
# print (rootProb)
# mu = -1*((rootProb[0]*Q[0,0])+(rootProb[1]*Q[1,1])+(rootProb[2]*Q[2,2])+(rootProb[3]*Q[3,3]))
# Q/=mu