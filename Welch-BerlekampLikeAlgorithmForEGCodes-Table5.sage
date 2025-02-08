def Encoding_Gabidulin(Message, Gabidulin_Support):
    f = S(Message.list())             #  The message polynomial 
    return vector(f.multi_point_evaluation(Gabidulin_Support))

def random_small_vec_gen(n,t): 
    B = matrix(Fqm.base_ring(),t,m,0)
    while B.rank() != t:
        B = matrix(Fqm.base_ring(),[vector(Fqm.random_element()) for i in range(t)])
    C = matrix(Fqm.base_ring(),n,m,0)
    while C.rank() != t: 
        C = matrix(Fqm.base_ring(),n,t,[Fqm.base_ring().random_element() for _ in range(n*t)]) * B  
    return vector(Fqm,[C[i] for i in range(n)])

def Interpolation_Polynomial(support_vector, recieved_word, dimension):
    gg = [support_vector[i] for i in range(dimension)]; yy = [recieved_word[i] for i in range(dimension)]
    Sum = 0
    for i in range(dimension): 
        temp = copy(gg); del temp[i]
        AP = S.minimal_vanishing_polynomial(temp)
        uu = AP(gg[i]).inverse_of_unit()
        Sum = Sum + yy[i] * uu * AP 
    return Sum

def swap_list_elements(lst, index1, index2):
    temp = lst[index1]
    lst[index1] = lst[index2]
    lst[index2] = temp
    return lst

def LRReconstructionAlgorithm(n, k, g, y):
    # Initialization Step
    gg = [g[i] for i in range(k)]
    N_0 = S.minimal_vanishing_polynomial(gg)
    W_0 = S.zero()
    N_1 = Interpolation_Polynomial(g, y, k)
    W_1 = S.one()
    
    u_0 = [N_0(g[i]) - W_0(y[i]) for i in range(n)]
    u_1 = [N_1(g[i]) - W_1(y[i]) for i in range(n)]
            
    # Interpolation Step
    for i in range(k, n):
        j = i
        if u_0[i] != 0 and u_1[i] == 0:
            for s in range(i, n):
                if u_1[s] == 0:
                    j = j + 1
                if u_1[s] != 0:
                    j = j 
            if j == n:
                return N_1, W_1
            else:
                #j = j
                u_0 = swap_list_elements(u_0, i, j)
                u_1 = swap_list_elements(u_1, i, j)
        
        # Updates of theta-polynomial, according discrepancies
        if u_1[i] != 0:
            N1 = (x - Frob(u_1[i]) * u_1[i].inverse_of_unit()) * N_1
            W1 = (x - Frob(u_1[i]) * u_1[i].inverse_of_unit()) * W_1
            N0 = N_0 - u_0[i] * u_1[i].inverse_of_unit() * N_1
            W0 = W_0 - u_0[i] * u_1[i].inverse_of_unit() * W_1
        
        if u_0[i] == 0 and u_1[i] == 0:
            N1 = x * N_1
            W1 = x * W_1
            N0 =  N_0
            W0 =  W_0
            
        N_0 =  N1 
        W_0 =  W1
        N_1 =  N0 
        W_1 =  W0
        
        # Updates of discrepancies
        u_0 = [N_0(g[i]) - W_0(y[i]) for i in range(n)]
        u_1 = [N_1(g[i]) - W_1(y[i]) for i in range(n)]
    return N_1, W_1


def WelchBerlekampDecoding(n, k, g, y):
    F = LRReconstructionAlgorithm(n, k, g, y)
    N_1 = F[0]; W_1 = F[1]
    ff, re = N_1.left_quo_rem(W_1)
    return ff

def test(totalltests):
    succ = 0
    failure = 0
    for npair in range(totalltests):
        e = random_small_vec_gen(n, r)
        y = Codeword + e
        ff = WelchBerlekampDecoding(n, k, g, y)
        try:
            if (vector(ff.padded_list(k))==Message) and (y - vector(ff.multi_point_evaluation(g)) == e): 
                succ += 1
            else:
                failure += 1
        except:
            print("Unexpected error", sys.exc_info()[0])
            
    print ("success/totalltests: %d/%d; success rate: %f" % (succ,totalltests,succ/totalltests))
    print ("failure/totalltests: %d/%d; failure rate: %f" % (failure,totalltests,failure/totalltests))

# Compute Theoretical and Simulated DFR by Theorem 3 for code parameters in Table 5
# increase m
#(q,m,n,t,k,r) = (2,31,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,32,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,33,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,34,41,34,9,16)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,35,41,35,9,16)   # DFR: 2**(-9)
# TheoreticalDFR = [0.0313, 0.0156, 0.0078, 0.0039, 0.0020]
# SimulatedDFR = [0.0150, 0.0076, 0.0038, 0.0017, 0.0008]

# increase t
#(q,m,n,t,k,r) = (2,35,41,30,9,16)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,35,41,31,9,16)   # DFR: 2**(-5)
#(q,m,n,t,k,r) = (2,35,41,32,9,16)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,35,41,33,9,16)   # DFR: 2**(-7)
#(q,m,n,t,k,r) = (2,35,41,34,9,16)   # DFR: 2**(-8)
# TheoreticalDFR = [0.0625, 0.0313, 0.0156, 0.0078, 0.0039]
# SimulatedDFR = [0.0310, 0.0150, 0.0084, 0.0036, 0.0018] 

# increase t < n < m; increase t
#(q,m,n,t,k,r) = (2,29,26,16,5,10)   # DFR: 2**(-2)
#(q,m,n,t,k,r) = (2,29,26,17,5,10)   # DFR: 2**(-4)
#(q,m,n,t,k,r) = (2,29,26,18,5,10)   # DFR: 2**(-6)
#(q,m,n,t,k,r) = (2,29,26,19,5,10)   # DFR: 2**(-8)
#(q,m,n,t,k,r) = (2,29,26,20,5,10)   # DFR: 2**(-10)
# TheoreticalDFR = [0.2500, 0.0625, 0.0156, 0.0039, 0.00098]
# SimulatedDFR = [0.1320, 0.0372, 0.0098, 0.0025, 0.00058] 

# Decoding up to the RGV bound (k > r)
# Hash-Sign 
#(q,m,n,t,k,r) = (2,30,37,30,23,7)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,38,30,22,8)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,39,30,21,9)   # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,30,40,30,20,10)  # DFR: 2**(1) 
#(q,m,n,t,k,r) = (2,30,41,30,19,11)  # DFR: 2**(1)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7120, 0.7087, 0.7100, 0.7105, 0.7141] 

# Decoding up to the RGV bound (k < r)
#(q,m,n,t,k,r) = (2,21,34,21,8,13)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,22,36,22,8,14)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,23,38,23,8,15)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,24,40,24,8,16)  # DFR: 2**(1)
#(q,m,n,t,k,r) = (2,25,42,25,8,17)  # DFR: 2**(1)
# TheoreticalDFR = [1, 1, 1, 1, 1]
# SimulatedDFR = [0.7124, 0.7120, 0.7118, 0.7120, 0.7125]

Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
# S = OrePolynomialRing(Fqm, Frob, 'x')
S.<x> = Fqm['x', Frob]

Message = random_vector(Fqm, k) 
g = random_small_vec_gen(n, min(m, n, t))
#g1 = random_small_vec_gen(t, t);  g2 = zero_vector(Fqm, n-t);  g = vector(g1.list() + g2.list())
Codeword = Encoding_Gabidulin(Message, g)


%time test(100000)
