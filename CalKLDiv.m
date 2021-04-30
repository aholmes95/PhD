function KLDiv = CalKLDiv(FP, G)
    %FP=[fp1, fp2, fp3] are the values of each household state respectively at
    %steady state (2,0), (1,1), (0,2) under the desired approximation method, 
    %while G=[g1, g2, g3] correspond to simulations from the Gillespie algorithm.
    
    %For discrete probability distributions P, Q, the KLDiv of Q (approximation) from
    %P (gillespie) is defined to be sum(Plog(P/Q)).
    
    Q = FP;
    P = G;
  
    KLDiv = sum(P.*log(P./Q));