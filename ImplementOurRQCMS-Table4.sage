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


# (q, m, n, t, k, r) = (2,83,79, 83, 7,36) 
(q, m, n, t, k, r) = (2, 53, 89, 53, 5, 1) 
(q, m, n, t, k, r) = (2, 5, 7, 5, 2, 2) 

Fqm.<a> = GF(q**m)
Frob = Fqm.frobenius_endomorphism()
# S = OrePolynomialRing(Fqm, Frob, 'x')
S.<x> = Fqm['x', Frob]

Message = random_vector(Fqm, k) 
g = random_small_vec_gen(n, min(m, n, t))
#g1 = random_small_vec_gen(t, t);  g2 = zero_vector(Fqm, n-t);  g = vector(g1.list() + g2.list())
Codeword = Encoding_Gabidulin(Message, g)
e = random_small_vec_gen(n,r)
y = Codeword + e

%time ff = WelchBerlekampDecoding(n, k, g, y)
    
print("Ture or False ? : ", y - vector(ff.multi_point_evaluation(g))== e and vector(ff.padded_list(k))==Message) 
