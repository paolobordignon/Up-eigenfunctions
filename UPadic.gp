

timeon = 1;
print("\n \tBefore starting select a value for ppp \n");
/*
EigenpV(n,ppp,i_eig,r) is a function to compute the eigenvalues of the Up operator 
and return the valuation of its zeroes wrt the parameter h=(Delta(px)/Delta(x))^(1/p+1). 
    n = precision of q expansion; 
    ppp = prime p; 
    i_eig = i-th eigenvector ordered by slope; 
    r = overconvergent parameter, r must be less then p/p+1 and 
        the denominator has to divide 12/(p-1)
*/

Eigenfunction(n,i_eig,r)={

[Mp, vecRoot] = Matrix_UP(n,1);

EigenpV(n,i_eig,r,Mp*(1+O(ppp^70)),vecRoot[i_eig]);

return(0);

}


EigenpV(n,i_eig,r,Mp,i_eigenvalue)={

    tott = getabstime();
    prec = n*5;
    k = 5; \\ 2^k is the power of the (Up-lambda Id) computed


    t=getabstime();

    slope_i = valuation(i_eigenvalue,ppp);
    \\V = ((Mp-lift(i_eigenvalue))/ppp^(vallambdai+1))^(-1);
    first_app_ev = lift(Mod(i_eigenvalue/ppp^slope_i,ppp));
    V = ((Mp-first_app_ev*ppp^slope_i)/ppp^(slope_i+1))^(-1);

    t=getabstime()-t;

    if(timeon!=0,
        print("____________________________________________");
        print("Time for computing inverse of Up matrix: "strtime(t));
        print("____________________________________________\n");
    );
    \\ for(i=1,l,V=(V*V));

    t = getabstime();

    for(l=1,k,V=V*V);

    t = getabstime()-t;

    if(timeon!=0,
        print("____________________________________________");
        print("Time for computing ",2^k," power of Up matrix: "strtime(t));
        print("____________________________________________");
        print("");
        print("");
    );
        eigenfunction_i = V[1,];
        eigenfunction_i = eigenfunction_i* ppp^(-6)/eigenfunction_i[1];

        precision_eigenfunction = valuation((eigenfunction_i*Mp-eigenfunction_i*i_eigenvalue)[1],ppp);

        print("The precision of the ",i_eig," eigenvector is: \t O(",ppp,"^",precision_eigenfunction,")\n");

    \\eigenfunction_i += vector(#eigenfunction_i,i_init_vec,O(ppp^precision_eigenfunction));

    eigenfunction_i = Clean_precision_vector(eigenfunction_i,precision_eigenfunction);

    NP_roots_ef = newtonpoly(Polrev(eigenfunction_i),ppp);

    for(i=1, #NP_roots_ef, NP_roots_ef[i]=-NP_roots_ef[i]*(ppp-1)/12+r);
    \\for(i=1, #NP_roots_ef, NP_roots_ef[i]=NP_roots_ef[i]-12*r/(ppp-1));
    
    NP_roots_ef = vecmultiplicities(NP_roots_ef);
    
    print("____________________________________________");
    for(i=1,#NP_roots_ef, 

        print("\n",NP_roots_ef[i][1],"\t w/ multiplicity \t",NP_roots_ef[i][2],"\n");

    );
    print("____________________________________________");

    tott = getabstime() - tott;

    print("\tTotal time: \t"strtime(tott));
    print("____________________________________________");

    return([slope_i,eigenfunction_i,precision_eigenfunction,NP_roots_ef]);
    
}

addhelp(EigenpV, "\n EigenpV(n,ppp,i_eig,r) is a function to compute the eigenvalues of the Up operator,\n return the valuation of its zeroes wrt the parameter h=(Delta(px)/Delta(x))^(1/p+1). \n\n\tn = precision of q expansion; \n\tppp = prime p; \n\ti_eig = i-th eigenvector ordered by slope; \n\tr = overconvergent parameter, \n\n Warning: the denominator of r has to divide 12/(p-1)");

Matrix_UP(n,r) = {

prec = n*10;


M2 = [48, 1;
    4096, 0];
M3 = [ 270 , 36 , 1 ;
    26244 , 729 , 0 ;
    531441, 0 , 0];
M5 = [ 1575, 1300 ,      315 ,       30 ,        1;
    162500 ,    39375 ,     3750 ,      125  ,       0;
    4921875 ,   468750 ,    15625 ,        0  ,       0;
    58593750 ,  1953125 ,        0 ,        0  ,       0;
    244140625 ,        0 ,        0 ,        0  ,       0];

M7 = [       4018 ,       8624 ,       5915 ,       1904 ,        322 ,         28 ,          1;
         422576 ,     289835 ,      93296  ,     15778      ,  1372      ,    49      ,     0;
       14201915  ,   4571504  ,    773122  ,     67228      ,  2401      ,     0      ,     0;
      224003696  ,  37882978  ,   3294172  ,    117649      ,     0      ,     0      ,     0;
     1856265922  , 161414428  ,   5764801  ,         0      ,     0      ,     0      ,     0;
     7909306972  , 282475249  ,         0  ,        0       ,   0        ,   0        ,   0;
    13841287201  ,         0  ,         0  ,        0       ,    0       ,    0       ,    0];

M13 = [   15145,         124852,         354536,         534820,         509366,         333580,         157118,          54340,          13832,           2548,            325,             26,              1;
        1623076,        4608968,        6952660,        6621758,       4336540 ,       2042534 ,        706420 ,        179816 ,         33124 ,          4225 ,           338 ,            13 ,             0;
       59916584,       90384580,       86082854,       56375020,       26552942,        9183460,        2337608,         430612,          54925,           4394,            169,              0,              0;
     1174999540,     1119077102,      732875260,      345188246,      119384980,       30388904,        5597956,         714025,          57122,           2197,              0,              0,              0;
    14548002326,     9527378380,     4487447198,     1552004740,      395055752,       72773428,        9282325,         742586,          28561,              0,              0,              0,              0;
    123855918940,    58336813574,    20176061620,     5135724776,      946054564,      120670225,        9653618,         371293,              0,              0,              0,              0,              0;
    758378576462,   262288801060,    66764422088,    12298709332,     1568712925,      125497034,        4826809,              0,              0,              0,              0,              0,              0;
    3409754413780,   867937487144,   159883221316,    20393268025,     1631461442,       62748517,              0,              0,              0,              0,              0,              0,              0;
    11283187332872,  2078481877108,   265112484325,    21208998746,      815730721,              0,              0,              0,              0,              0,              0,              0,              0;
    27020264402404,  3446462296225,   275716983698,    10604499373,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    44804009850925,  3584320788074,   137858491849,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    46596170244962,  1792160394037,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    23298085122481,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0];


Ip2rec = 1-[x,x^2]*M2*mattranspose([y,y^2]);
Ip3rec = 1-[x,x^2,x^3]*M3*mattranspose([y,y^2,y^3]);
Ip5rec = 1-[x,x^2,x^3,x^4,x^5]*M5*mattranspose([y,y^2,y^3,y^4,y^5]);
Ip7rec = 1-[x,x^2,x^3,x^4,x^5,x^6,x^7]*M7*mattranspose([y,y^2,y^3,y^4,y^5,y^6,y^7]);
Ip13rec = 1-[x,x^2,x^3,x^4,x^5,x^6,x^7,x^8,x^9,x^10,x^11,x^12,x^13]*M13*mattranspose([y,y^2,y^3,y^4,y^5,y^6,y^7,y^8,y^9,y^10,y^11,y^12,y^13]);

Mprec = matrix(ppp);

if(ppp == 2,
    Ip = Ip2rec;
    Mprec = M2;
);
if(ppp == 3,
    Ip = Ip3rec;
    Mprec = M3;
);
if(ppp == 5,
    Ip = Ip5rec;
    Mprec = M5;
);
if(ppp == 7,
    Ip = Ip7rec;
    Mprec = M7;
);
if(ppp == 13,
    Ip = Ip13rec;
    Mprec = M13;
);

t=getabstime();

f =Pol( deriv(Ip,y)/Ip*y/ppp+O(x^(ppp+1)));
Mp = matrix(n,n,x,y,O(ppp^prec));

t=getabstime()-t;

if(timeon!=0,
    print("____________________________________________");
    print("Time for computing power series of f: ",strtime(t));
    print("____________________________________________");
    print("");
    print("");
);

c = ppp^(12*r/(ppp-1));

\\ a base for X_0(p)[r] is given by cf, (cf)^2, (cf)^3, ...

t=getabstime();

Vcoef = -1*Vecrev(f);
for(i=1,ppp,Vcoef[i]=Pol(Vcoef[i]+O(y^(ppp+1))));

for(j=2,ppp+1,
    for(i=2,min(#Vec(Vcoef[j]),ppp+1), 
        Mp[i-1,j-1]= Vecrev(Vcoef[j])[i];
        );
    );   

for(i=1, n,
    for(j=1,n,
        if(i > ppp || j > ppp, 
            for(a=1,min(ppp,j-1),
            for(b=1,min(ppp,i-1),
                Mp[i,j]=Mp[i,j]+Mprec[a,b]*Mp[i-b,j-a];
            );
        );
        );
    );
);

\\ complexity p^2*n^2

for(i=1, n,
    for(j=1,n,
        Mp[i,j]=Mp[i,j]*c^(i-j);
    );
);

t=getabstime()-t;

if(timeon!=0,
    print("____________________________________________");
    print("Time for computing Up matrix: ",strtime(t));
    print("____________________________________________");
    print("");
    print("");
);

/*
for(i=1,n,for(j=1,n, 
    if(Mp[i,j]!=0,vMp[i,j]=valuation(Mp[i,j],ppp));
    ));
*/

    Mp = Mp*(1+O(ppp^prec));

    t=getabstime();
    char_pol_Up = charpoly(Mp);
    t=getabstime()-t;
    
    if(timeon!=0,
        print("");
        print("____________________________________________");
        print("Time for computing characteristic polynomial: ",strtime(t));
        print("____________________________________________");
        print("");
        print("");
    );

    NP_Up = newtonpoly(char_pol_Up,ppp); 
    vecRoot = Vecrev(polrootspadic(char_pol_Up,ppp,prec));

    print("\n","\t",#vecRoot," eigenvalues computed\n");

    return([Mp,vecRoot]);

}

vecmultiplicities(v) = {

    vecEntries = List();
    firstinstance = 1;

    for(i=2,#v,
        if(v[i]!=v[i-1],
        listput(vecEntries,[v[firstinstance],i-firstinstance]);
        firstinstance = i;
        );
    );
    listput(vecEntries,[v[#v],#v-firstinstance+1]);
    return(Vec(vecEntries));
    
}

etaqexpansion(phivect,r)={

    c = ppp^(12*r/(ppp-1));
    sprec = 100;
    etaq = 1+O(q^sprec);

    for(i=1,sprec, etaq= etaq*(1-q^i));

    etaqp = subst(etaq,q,q^ppp);
    hmodul =c*(etaqp/etaq)^(24/(ppp-1))*q;

    phiqexp = O(q^sprec);
    
    for(i = 1, min(sprec,#phivect), phiqexp+=hmodul^i*phivect[i]);

    return(phiqexp/Vec(phiqexp)[1]);
}


/*

See [Smi07] "Bounding slopes of p-adic modular forms" or [Loe07] "Spectral expansion of overconvergent modular functions"
We want to compute a polynomial F_p(x,y) such that F_p(V(h),h)=0. 
This should play the role of a planar model for X_0(p^2)
*/


Fppol()={

f = Recursive_Matrix(20)[1];

Fp = subst(f,y,1/y)*y^ppp;

Fp_matrix = matrix(ppp+1);

for(i=1, min(#Vecrev(Fp),ppp+1),
    for(j=1,min(#Vecrev(Vecrev(Fp)[i]),ppp+1),
                Fp_matrix[i,j]=Vecrev(Vecrev(Fp)[i])[j];
    );
);

return(Fp);

}

ZeroeslogE0(depth)={

    Fp = Fppol();

    print("____________________________________________");
    print("Valuation of zeroes of E_0^(",ppp,")");
    print("____________________________________________");
    print("");
    print("");
    for(i=-depth,depth, 
        Gp = subst(Fp,x,y^ppp*ppp^i);
        inford = valuation(Gp,y);
        print("____________________________________________");
    print("The zeroes of F_",ppp,"(x,x^",ppp,"*",ppp,"^",i,")");
    print("");
        print(newtonpoly(Gp/y^inford,ppp)*(ppp-1)/(-12));
        print("");
    print("____________________________________________");
    print("");
    print("");
    );

}

StoreUpeigenfunction(L_store,i1,i2,n,padic_prec)={

    timeon = 0;

    [Mp, vecRoot] = Matrix_UP(n,1);

    for(i=i1,i2,
    print("____________________________________________");
    print("\nComputing the ",i,"-th eigenfunction...\n");
    print("____________________________________________");

    i_eigenvalue = vecRoot[i];

    info_i_eig = EigenpV(n,i,1,Mp*(1+O(ppp^padic_prec)),i_eigenvalue); 
    slope_i = info_i_eig[1];
    eigenfunction_i = info_i_eig[2];
    precision_eigenfunction = info_i_eig[3];

    listput(L_store,[Str(i,"-th ",ppp,"-adic eigenfunction for U_",ppp," operator with slope ",slope_i, " (radius of convergence 1)"),info_i_eig],i);

    );

    print("\n The returned list contains in the second argument \n \t 1. slope,\n\t 2. eigenfunction with respect to basis h,\n\t 3. precision,\n\t 4. slopes of Newton Polygon wrt basis h");

    return(L_store);

}

Clean_precision_vector(vec_V,prec_V)={

    len_V = #vec_V;
    \\vec_V += vector(len_V,i_init_vec,O(ppp^prec_V));
    vec_V = vec_V*(1+O(ppp^prec_V));

    for(i=0,len_V-1,
        if(lift(vec_V[len_V-i])==0, vec_V[len_V-i]=0, i = len_V);
    );

    return(vec_V);

}


Up_to_Fp(n)={

[Ip,Mp_rec] = Recursive_Matrix(n);


f = Pol( deriv(Ip,y)/Ip*y/ppp+O(x^(n+1)));

Up_f = matrix(n);
Up_val = matrix(n);

Ip_matrix = matrix(ppp+1);
Ip_matrix_val = matrix(ppp+1);

f = f*(-1);

for(i=1, min(#Vecrev(f)-1,n),
    for(j=1,min(#Vecrev(Vecrev(f)[i+1])-1,n),
                Up_f[i,j]=Vecrev(Vecrev(f)[i+1])[j+1];
    );
);

for(i=1,n,
    for(j=1,n,
        Up_val[i,j] = valuation(Up_f[i,j],ppp);
        if(Up_f[i,j]==0,Up_val[i,j]=x);
    );
);

for(i=1, min(#Vecrev(Ip)-1,ppp+1),
    for(j=1,min(#Vecrev(Vecrev(Ip)[i+1])-1,ppp+1),
                Ip_matrix[i,j]=Vecrev(Vecrev(Ip)[i+1])[j+1];
    );
);

for(i=1,ppp+1,
    for(j=1,ppp+1,
        Ip_matrix_val[i,j] = valuation(Ip_matrix[i,j],ppp);
        if(Ip_matrix[i,j]==0,Ip_matrix_val[i,j]=x);
    );
);

\\ complexity p^2*n^2


/*
for(i=1, n,
    for(j=1,n,
        Mp[i,j]=Mp[i,j]*ppp^(i-j);
    );
);
*/

g = 0;

for(i=1, n,
    for(j=1,n,
        g += Up_f[i,j]*ppp/j*x^i*y^j;
    );
);


bound_matrix = matrix(n);

for(i=1, n,
    for(j=1,n,
        bound_matrix[i,j]=floor(12*(ppp*i-j)/((ppp+1)*(ppp-1)))-1;
    );
);

return(Up_f);

}

Recursive_Matrix(n)={

M2 = [48, 1;
    4096, 0];
M3 = [ 270 , 36 , 1 ;
    26244 , 729 , 0 ;
    531441, 0 , 0];
M5 = [ 1575, 1300 ,      315 ,       30 ,        1;
    162500 ,    39375 ,     3750 ,      125  ,       0;
    4921875 ,   468750 ,    15625 ,        0  ,       0;
    58593750 ,  1953125 ,        0 ,        0  ,       0;
    244140625 ,        0 ,        0 ,        0  ,       0];

M7 = [       4018 ,       8624 ,       5915 ,       1904 ,        322 ,         28 ,          1;
         422576 ,     289835 ,      93296  ,     15778      ,  1372      ,    49      ,     0;
       14201915  ,   4571504  ,    773122  ,     67228      ,  2401      ,     0      ,     0;
      224003696  ,  37882978  ,   3294172  ,    117649      ,     0      ,     0      ,     0;
     1856265922  , 161414428  ,   5764801  ,         0      ,     0      ,     0      ,     0;
     7909306972  , 282475249  ,         0  ,        0       ,   0        ,   0        ,   0;
    13841287201  ,         0  ,         0  ,        0       ,    0       ,    0       ,    0];

M13 = [   15145,         124852,         354536,         534820,         509366,         333580,         157118,          54340,          13832,           2548,            325,             26,              1;
        1623076,        4608968,        6952660,        6621758,       4336540 ,       2042534 ,        706420 ,        179816 ,         33124 ,          4225 ,           338 ,            13 ,             0;
       59916584,       90384580,       86082854,       56375020,       26552942,        9183460,        2337608,         430612,          54925,           4394,            169,              0,              0;
     1174999540,     1119077102,      732875260,      345188246,      119384980,       30388904,        5597956,         714025,          57122,           2197,              0,              0,              0;
    14548002326,     9527378380,     4487447198,     1552004740,      395055752,       72773428,        9282325,         742586,          28561,              0,              0,              0,              0;
    123855918940,    58336813574,    20176061620,     5135724776,      946054564,      120670225,        9653618,         371293,              0,              0,              0,              0,              0;
    758378576462,   262288801060,    66764422088,    12298709332,     1568712925,      125497034,        4826809,              0,              0,              0,              0,              0,              0;
    3409754413780,   867937487144,   159883221316,    20393268025,     1631461442,       62748517,              0,              0,              0,              0,              0,              0,              0;
    11283187332872,  2078481877108,   265112484325,    21208998746,      815730721,              0,              0,              0,              0,              0,              0,              0,              0;
    27020264402404,  3446462296225,   275716983698,    10604499373,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    44804009850925,  3584320788074,   137858491849,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    46596170244962,  1792160394037,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0;
    23298085122481,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0];


Ip2rec = 1-[x,x^2]*M2*mattranspose([y,y^2]);
Ip3rec = 1-[x,x^2,x^3]*M3*mattranspose([y,y^2,y^3]);
Ip5rec = 1-[x,x^2,x^3,x^4,x^5]*M5*mattranspose([y,y^2,y^3,y^4,y^5]);
Ip7rec = 1-[x,x^2,x^3,x^4,x^5,x^6,x^7]*M7*mattranspose([y,y^2,y^3,y^4,y^5,y^6,y^7]);
Ip13rec = 1-[x,x^2,x^3,x^4,x^5,x^6,x^7,x^8,x^9,x^10,x^11,x^12,x^13]*M13*mattranspose([y,y^2,y^3,y^4,y^5,y^6,y^7,y^8,y^9,y^10,y^11,y^12,y^13]);

Mprec = matrix(ppp);

if(ppp == 2,
    Ip = Ip2rec;
    Mprec = M2;
);
if(ppp == 3,
    Ip = Ip3rec;
    Mprec = M3;
);
if(ppp == 5,
    Ip = Ip5rec;
    Mprec = M5;
);
if(ppp == 7,
    Ip = Ip7rec;
    Mprec = M7;
);
if(ppp == 13,
    Ip = Ip13rec;
    Mprec = M13;
);


f = Pol( deriv(Ip,y)/Ip*y/ppp+O(x^(n+1)));

return([Ip,Mprec]);

}

sym_Up(n)={

    sym = matrix(n,n,i,j,if(i==j,i));

    return(sym);
    
}