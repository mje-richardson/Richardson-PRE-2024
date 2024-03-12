#########################################################################
#########################################################################
#########################################################################
# Contains the parameters and functions that are common to the various figures
#
#########################################################################
#########################################################################
# --- COMMON PARAMETERS ---
#########################################################################
#########################################################################
#
#########################################################################
#########################################################################
# --- GENERAL FUNCTIONS ---
#########################################################################
#########################################################################
#               
# v=MyfRoots(v1,vT,dT,n)       
# r0,amp,phi=MyExtractFourier(w,t,Rs1,rst)
#
#########################################################################
#########################################################################
# --- THEORY FUNCTIONS ---
#########################################################################
#########################################################################
#
#########################################################################
# ThIn0 functions 
#########################################################################
# r,v,P,Je,Ji=LIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau)
# r,v,P,Je,Ji=EIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT)  
# r,v,P,Je,Ji=LIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,ee,ei)
# r,v,P,Je,Ji=EIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,ee,ei)
#
#########################################################################
# ThIn1 functions
#########################################################################
# r,re,ri=LIFShotCurrentThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,fkHz)
# r,re,ri=EIFShotCurrentThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,fkHz)
# r,re,ri=LIFShotConductThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,fkHz,ee,ei)
# r,re,ri=EIFShotConductThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,fkHz,ee,ei)
#
#########################################################################
#########################################################################
# --- SIMULATION FUNCTIONS ---
#########################################################################
#########################################################################
#
# v,vend,r,rsum=EIFShotCurrentSimEI(tau,dT,vT,vth,vre,vlb,t,Ret,Rit,ae,ai,v1)
# v,vend,r,rsum=EIFShotCurrentSimEI(tau,dT,vT,vth,vre,vlb,t,Ret,Rit,ae,ai,v1)
#
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# Common parameters across the figures
#########################################################################
#########################################################################
#########################################################################
#########################################################################

########################################################################
# sundry parameters
#########################################################################
K=1000
r2d=360/(2*pi)

#########################################################################
# shared IF parameters
#########################################################################
tau=20
vre=5
ee=60
ei=-10
vlbCond=ei+0.1;
vlbCurr=-15.0;

#########################################################################
# EIF neuron parameters
#########################################################################
vT=10
dT=1
vthEIF=20

#########################################################################
# LIF neuron parameters
#########################################################################
vthLIF=10

#########################################################################
# synaptic amplitudes and reversal potentials
#########################################################################
ae=1.5
ai=-0.75
sig=5

#########################################################################
# numerical parameters
#########################################################################
dv=0.0005 # for ThIn method
dt=0.01 # for simulations

println("Neuron params:\ttau=$(tau)ms, ee=$(ee)mV, ei=$(ei)mV vT=$(vT)mV, dT=$(dT)mV");
println("Thresh params:\tvthEIF=$(vthEIF)mV vthLIF=$(vthLIF)")
println("LowerB params:\tvlbCond=$(vlbCond)mV vlbCurr=$(vlbCurr)")
println("Syn params:\tae,ai=$ae, $(ai)mV and sig=$(sig)mV")
println("Num params:\tThIn dv=$(dv)mV & sims dt=$(dt)ms")

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# GENERAL FUNCTIONS
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# Halley's method for lower and upper fixed points. 
# Needed for LambertW calculate of EIF fixed points
#########################################################################
#########################################################################
function MyfRoots(v1,vT,dT,n)
   
    v=zeros(n)
    v[1]=v1
    for k=1:n-1
        expk=exp((v[k]-vT)/dT)
        f0=dT*expk-v[k]
        f1=expk-1
        f2=(1/dT)*expk
        v[k+1]=v[k]-2*f0*f1/(2*f1^2-f0*f2)
    end

    return v[end]
    
end

#########################################################################
#########################################################################
# Extract fourier amplitudes and phase
#########################################################################
#########################################################################
function MyExtractFourier(w,t,Rs1,rst)
    
    nw=length(w)
    r0=mean(rst) 
    amp,phi=zeros(nw),zeros(nw)
    
    for k=1:nw
        X1=(2/T)*sum(dt*rst.*cos.(w[k]*t))/Rs1 
        X2=(2/T)*sum(dt*rst.*sin.(w[k]*t))/Rs1
        amp[k]=sqrt(X1^2 + X2^2) 
        phi[k]=angle(X1 - 1im*X2)
    end
    
    return r0,amp,phi

end

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# STEADY-STATE FUNCTIONS
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# LIFshotEIcurrent0ThIn 
# Solution using Threshold Integration method for the 
# LIF driven by current-based shot noise
#########################################################################
#########################################################################
function LIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau)
    
    v=collect(vlb+dv/2:dv:vth-dv/2)    
    n=length(v)
    F=-v/tau
    
    # find points above and below the origin
    n0p=findall(v .>0)[1]
    n0m=n0p-1
    
    # initialise the arrays
    (p,je,ji)=(zeros(n) for i=1:3)
    
    ###############################################################
    # integrate over the range      vlb >-----> 0 ....vre..... vth
    ###############################################################
    
    je[1]=0 
    ji[1]=-1
    p[1]=-(je[1]+ji[1])./F[1]
    
    for k=1:1:n0m-1
        je[k+1]=je[k]+dv*(Re*p[k]-je[k]/ae)
        ji[k+1]=ji[k]+dv*(Ri*p[k]-ji[k]/ai)   
        p[k+1]=-(je[k+1]+ji[k+1])./F[k+1]
    end
    
    ###############################################################
    # integrate over the range      vlb ........ 0 <---vre---< vth
    ###############################################################
    
    je[n]=1 
    ji[n]=0 
    p[n]=((v[n].>vre)-je[n]-ji[n])./F[n];
    
    for k=n:-1:n0p+1
        je[k-1]=je[k]-dv*(Re*p[k]-je[k]/ae)
        ji[k-1]=ji[k]-dv*(Ri*p[k]-ji[k]/ai)
        p[k-1]=((v[k-1].>vre)-je[k-1]-ji[k-1])./F[k-1]
    end
    
    # match v<0 and v>0 solns at v=0 using the excitatory flux
    a=je[n0p]/je[n0m]
    for f in [p,je,ji] f[1:n0m]=a*f[1:n0m] end
    
    # normalizing the distribution
    r=1/sum(p*dv);  
    P=r*p;
    Je=r*je;
    Ji=r*ji;

    return r,v,P,Je,Ji    

end

#########################################################################
#########################################################################
# EIFshotEIcurrent0ThIn 
# Solution using Threshold Integration method for the 
# EIF driven by current-based shot noise
#########################################################################
#########################################################################
function EIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT)
    
    # find the stable and unstable fixed points
    vs,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)
    
    # set things up
    v=collect(vlb-dv/2:dv:vth-dv/2)
    n=length(v)
    F=(dT*exp.((v .-vT)/dT) .-v)/tau
    
    # the reset
    krep=findall(v.>0)[1]
    krem=krep-1
    
    # points just above and below stable fixed point
    nsp=findall(v.>vs)[1]
    nsm=nsp-1
    
    # points above and below unstable fixed point
    nup=findall(v.>vu)[1]
    num=nup-1
    dFdvu=(exp((vu-vT)/dT)-1)/tau
    dup=v[nup]-vu
    dum=vu-v[num]
    
    # intitialise the arrays
    (p,pA,pB, je,jeA,jeB, ji,jiA,jiB)=(zeros(n) for i=1:9)
 
    ################################################################
    # integrate over the range   vlb >---> vs ...vre ... vu ... vth
    ################################################################
    je[1]=0; ji[1]=-1; p[1]=-(je[1]+ji[1])./F[1];
    for k=1:1:nsm-1
        je[k+1]=je[k]+dv*(Re*p[k]-je[k]/ae)
        ji[k+1]=ji[k]+dv*(Ri*p[k]-ji[k]/ai)
        p[k+1]=-(je[k+1]+ji[k+1])./F[k+1]
    end
    
    ################################################################
    # integrate over the range   vlb ... vs ...vre ... vu >---> vth
    ################################################################

    ##################################################################
    # we need to integrate two solutions from vu
    #
    # solution A has conditions jA=1 jeA=1 and jiA=0
    # solution B has conditions jB=0 jeB=1 and jiB=-1
    # and use the condition that ji(vth)=0 to find the superposition
    #
    # good estimates of the initial conditions just above the 
    # unstable fixed point are required for stability.
    ##################################################################
    
    jeAu=1; jiAu=0;  pAu=(jeAu/ae+jiAu/ai)/(dFdvu+Re+Ri);
    jeBu=1; jiBu=-1; pBu=(jeBu/ae+jiBu/ai)/(dFdvu+Re+Ri);
    
    # now improve the initial conditions - good for stability
    jeA[nup]=1
    jiA[nup]=0
    pA[nup]=(1-jeA[nup]-jiA[nup])./F[nup]
    
    jeB[nup]=1 
    jiB[nup]=-1
    pB[nup]=-(jeB[nup]+jiB[nup])./F[nup]
    
    for k=nup:n-1
        jeA[k+1]=jeA[k]+dv*(Re*pA[k]-jeA[k]/ae)
        jiA[k+1]=jiA[k]+dv*(Ri*pA[k]-jiA[k]/ai)
        pA[k+1]=(1-jeA[k+1]-jiA[k+1])./F[k+1]
        
        jeB[k+1]=jeB[k]+dv*(Re*pB[k]-jeB[k]/ae)
        jiB[k+1]=jiB[k]+dv*(Ri*pB[k]-jiB[k]/ai)
        pB[k+1]=-(jeB[k+1]+jiB[k+1])./F[k+1]
    end
    
    # this is the factor for the B solution
    b=-jiA[n]/jiB[n]
    
    # this is a first estimate and is improved below
    je[nup:n]=jeA[nup:n]+b*jeB[nup:n] 
    ji[nup:n]=jiA[nup:n]+b*jiB[nup:n] 
    p[nup:n]=pA[nup:n]+  b*pB[nup:n] 

    ###################################################################
    # integrate over the range   vlb ... vs <---< vre <---< vu ... vth
    ###################################################################
    
    jeu=jeAu+b*jeBu; jiu=jiAu+b*jiBu; pu=(jeu/ae+jiu/ai)/(dFdvu+Re+Ri)
    
    # now improve the initial conditions -  good for stability
    je[num]=jeu-dum*(Re*pu-jeu/ae) 
    ji[num]=jiu-dum*(Ri*pu-jiu/ai)
    p[num]=(1-je[num]-ji[num])./F[num]
    
    for k=num:-1:nsp+1
        je[k-1]=je[k]-dv*(Re*p[k]-je[k]/ae)
        ji[k-1]=ji[k]-dv*(Ri*p[k]-ji[k]/ai)
        p[k-1]=((v[k-1].>vre)-je[k-1]-ji[k-1])./F[k-1]
    end
  
    ##############################################################
    # match v<vs and v>vs solns at v=vs using the excitatory flux
    ##############################################################
    
    a=je[nsp]/je[nsm]
    for f in [p,je,ji] f[1:nsm]=a*f[1:nsm] end
    
    r=1/sum(p*dv)
    P=r*p
    Je=r*je;
    Ji=r*ji;

    return r,v,P,Je,Ji
    
end
    
#########################################################################
#########################################################################
# LIFshotEIconduct0ThIn 
# Solution using Threshold Integration method for the 
# LIF driven by current-based shot noise
#########################################################################
#########################################################################
function LIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,ee,ei)
    
    v=collect(vlb+dv/2:dv:vth-dv/2)
    n=length(v)
    F=-v/tau
    
    # find points above and below the origin
    n0p=findall(v.>0)[1]
    n0m=n0p-1
    
    # convert PSPs at v=0 to conductance parameters
    be=ae/ee
    bi=ai/ei
    betae=1/be-1
    betai=1/bi-1
    
    # intitialise the arrays
    (p,je,ji)=(zeros(n) for i=1:3)
    
    #########################################################################
    # integrate over the range      vlb >-----> 0 ....vre..... vth
    #########################################################################
    je[1]=0 
    ji[1]=-1 
    p[1]=-(je[1]+ji[1])./F[1]
    
    for k=1:1:n0m-1
        je[k+1]=je[k]+dv*(Re*p[k]-betae*je[k]/(ee-v[k]));
        ji[k+1]=ji[k]+dv*(Ri*p[k]-betai*ji[k]/(ei-v[k]))
        p[k+1]=-(je[k+1]+ji[k+1])./F[k+1];
    end
    
    #########################################################################
    # integrate over the range      vlb ........ 0 <---vre---< vth
    #########################################################################
    je[n]=1 
    ji[n]=0 
    p[n]=((v[n]>vre)-je[n]-ji[n])./F[n]
    
    for k=n:-1:n0p+1
        je[k-1]=je[k]-dv*(Re*p[k]-betae*je[k]/(ee-v[k]))
        ji[k-1]=ji[k]-dv*(Ri*p[k]-betai*ji[k]/(ei-v[k]))
        p[k-1]=((v[k-1]>vre)-je[k-1]-ji[k-1])./F[k-1]
    end
    
    #########################################################################
    # match v<0 and v>0 solns at v=0 using the excitatory flux
    #########################################################################
    a=je[n0p]/je[n0m]
    for f in [p,je,ji] f[1:n0m]=a*f[1:n0m] end
    
    # normalizing the distribution
    r=1/sum(p*dv)  
    P=r*p
    Je=r*je;
    Ji=r*ji;
    
    return r,v,P,Je,Ji

end                

#########################################################################
#########################################################################
# EIFshotEIcurrent0ThIn 
# Solution using Threshold Integration method for the 
# EIF driven by current-based shot noise
#########################################################################
#########################################################################
function EIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,ee,ei)

    # find the stable and unstable fixed points
    vs,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)
    
    #########################################################################
    # set things up
    #########################################################################
    v=collect(vlb-dv/2:dv:vth-dv/2)
    n=length(v)
    F=(dT*exp.((v .-vT)/dT)-v)/tau
    
    # synaptic parameters
    be=ae/ee
    bi=ai/ei
    betae=1/be-1
    betai=1/bi-1

    # the reset
    krep=findall(v.>0)[1]
    krem=krep-1

    # points just above and below stable fixed point
    nsp=findall(v.>vs)[1]
    nsm=nsp-1
    
    # points above and below unstable fixed point
    nup=findall(v.>vu)[1]
    num=nup-1
    dFdvu=(exp((vu-vT)/dT)-1)/tau
    dup=v[nup]-vu
    dum=vu-v[num]
    
    # intitialise the arrays
    (p,pA,pB, je,jeA,jeB, ji,jiA,jiB)=(zeros(n) for i=1:9)
    
    ########################################################################
    # integrate over the range   vlb >---> vs ...vre ... vu ... vth
    ########################################################################
    je[1]=0; ji[1]=-1; p[1]=-(je[1]+ji[1])./F[1]
    for k=1:1:nsm-1
        je[k+1]=je[k]+dv*(Re*p[k]-betae*je[k]/(ee-v[k]))
        ji[k+1]=ji[k]+dv*(Ri*p[k]-betai*ji[k]/(ei-v[k]))
        p[k+1]=-(je[k+1]+ji[k+1])./F[k+1]
    end
    
    ########################################################################
    # integrate over the range   vlb ... vs ...vre ... vu >---> vth
    ########################################################################
    # we need to integrate two solutions from vu
    #
    # solution A has conditions jA=1 jeA=1 and jiA=0
    # solution B has conditions jB=0 jeB=1 and jiB=-1
    # and use the condition that ji(vth)=0 to find the superposition
    #
    # good estimates of the initial conditions just above the 
    # unstable fixed point are required for stability.
    ########################################################################
    jeAu=1; jiAu=0;  pAu=(betae*jeAu/(ee-vu)+betai*jiAu/(ei-vu))/(dFdvu+Re+Ri);
    jeBu=1; jiBu=-1; pBu=(betae*jeBu/(ee-vu)+betai*jiBu/(ei-vu))/(dFdvu+Re+Ri);
    
    # now improve the initial conditions - good for stability
    jeA[nup]=1
    jiA[nup]=0
    pA[nup]=(1-jeA[nup]-jiA[nup])./F[nup];

    jeB[nup]=1
    jiB[nup]=-1
    pB[nup]=-(jeB[nup]+jiB[nup])./F[nup];
    
    for k=nup:n-1

        jeA[k+1]=jeA[k]+dv*(Re*pA[k]-betae*jeA[k]/(ee-v[k]))
        jiA[k+1]=jiA[k]+dv*(Ri*pA[k]-betai*jiA[k]/(ei-v[k]))
        pA[k+1]=(1-jeA[k+1]-jiA[k+1])./F[k+1]
        
        jeB[k+1]=jeB[k]+dv*(Re*pB[k]-betae*jeB[k]/(ee-v[k]))
        jiB[k+1]=jiB[k]+dv*(Ri*pB[k]-betai*jiB[k]/(ei-v[k]))
        pB[k+1]=-(jeB[k+1]+jiB[k+1])./F[k+1]
    
    end
    
    # this is the factor for the B solution
    b=-jiA[n]/jiB[n]
    
    # this is a first estimate and is improved below
    je[nup:n]=jeA[nup:n]+b*jeB[nup:n] 
    ji[nup:n]=jiA[nup:n]+b*jiB[nup:n] 
    p[nup:n]=pA[nup:n]+  b*pB[nup:n] 
    
    #########################################################################
    # integrate over the range   vlb ... vs <---< vre <---< vu ... vth
    #########################################################################
    jeu=jeAu+b*jeBu; jiu=jiAu+b*jiBu; 
    pu=(betae*jeu/(ee-vu)+betai*jiu/(ei-vu))/(dFdvu+Re+Ri)
    
    # now improve the initial conditions -  good for stability
    je[num]=jeu
    ji[num]=jiu
    p[num]=pu
    
    for k=num:-1:nsp+1
    
    je[k-1]=je[k]-dv*(Re*p[k]-betae*je[k]/(ee-v[k]))
    ji[k-1]=ji[k]-dv*(Ri*p[k]-betai*ji[k]/(ei-v[k]))   
    p[k-1]=((v[k-1]>vre)-je[k-1]-ji[k-1])./F[k-1]
    
    end

    #########################################################################
    # match v<vs and v>vs solns at v=vs using the excitatory flux
    #########################################################################
    a=je[nsp]/je[nsm]
    for f in [p,je,ji] f[1:nsm]=a*f[1:nsm] end
    
    r=1/sum(p*dv)
    P=r*p
    Je=r*je;
    Ji=r*ji;
    
    return r,v,P,Je,Ji

end

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# MODULATIONS FUNCTIONS
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# ThIn for modulated shot LIF with E and I 
#########################################################################
#########################################################################
function LIFShotCurrentThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,fkHz)
    
    # get useful steady-state parameters
    r,v,P,~,~=LIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau)
    
    v=collect(vlb+dv/2:dv:vth-dv/2)
    n=length(v)
    F=-v/tau
    
    # find points above and below the origin
    n0p=findall(v.>0)[1]
    n0m=n0p-1
    
    # the reset
    kre=minimum(findall(v.>vre))
    kreA=zeros(n)
    kreA[kre]=1
    
    # NOW THE MODULATION
    w=2*pi*fkHz;
    
    X=[dv*(1.0./F[1:n-1] .+ 1.0./F[2:n])/2;0]
    
    Xw=exp.(-1im*w*X); Yw=1.0 ./Xw;
    Xe=exp.(-Re*X);    Ye=1.0 ./Xe;
    Xi=exp.(-Ri*X);    Yi=1.0 ./Xi;

    # intitialise the arrays for the various solutions
    (sj,sp,sje,sji)=(zeros(Complex{Float64},n) for i=1:4)
    (rj,rp,rje,rji)=(zeros(Complex{Float64},n) for i=1:4)
    (ej,ep,eje,eji)=(zeros(Complex{Float64},n) for i=1:4)
    (ij,ip,ije,iji)=(zeros(Complex{Float64},n) for i=1:4)

    ##############################################################
    # integrate over the range      vlb >-----> 0 ....vre.... vth
    ##############################################################
    sj[1]=0; sje[1]=0; sji[1]=-1; # the boundary condition at vlb
    ej[1]=0; eje[1]=0; eji[1]=0;  
    ij[1]=0; ije[1]=0; iji[1]=0; 
    
    for k=1:1:n0m-1

        sj[k+1]=sj[k]*Xw[k]+(sje[k]+sji[k])*(1-Xw[k]);
        sje[k+1]=exp(-dv/ae)*(sje[k]*Xe[k]+(sj[k]-sji[k])*(1-Xe[k]));
        sji[k+1]=exp(-dv/ai)*(sji[k]*Xi[k]+(sj[k]-sje[k])*(1-Xi[k])); 

        ej[k+1]=ej[k]*Xw[k]+(eje[k]+eji[k])*(1-Xw[k]);
        eje[k+1]=exp(-dv/ae)*(eje[k]*Xe[k]+(ej[k]-eji[k])*(1-Xe[k]))+P[k]*dv;
        eji[k+1]=exp(-dv/ai)*(eji[k]*Xi[k]+(ej[k]-eje[k])*(1-Xi[k])); 

        ij[k+1]=ij[k]*Xw[k]+(ije[k]+iji[k])*(1-Xw[k]);
        ije[k+1]=exp(-dv/ae)*(ije[k]*Xe[k]+(ij[k]-iji[k])*(1-Xe[k]));
        iji[k+1]=exp(-dv/ai)*(iji[k]*Xi[k]+(ij[k]-ije[k])*(1-Xi[k]))+P[k]*dv
    
    end

    ###############################################################
    # integrate over the range      vlb ........ 0 <---vre---< vth
    ###############################################################
    rj[n]=1; rje[n]=1; rji[n]=0 # the boundary condition at vth
    ej[n]=0; eje[n]=0; eji[n]=0 
    ij[n]=0; ije[n]=0; iji[n]=0 
    
    for k=n:-1:n0p+1

        rj[k-1]=rj[k]*Yw[k]+(rje[k]+rji[k])*(1-Yw[k]) - kreA[k]
        rje[k-1]=exp(dv/ae)*rje[k]*Ye[k]+(rj[k]-rji[k])*(1-Ye[k])
        rji[k-1]=exp(dv/ai)*rji[k]*Yi[k]+(rj[k]-rje[k])*(1-Yi[k])
        
        ej[k-1]=ej[k]*Yw[k]+(eje[k]+eji[k])*(1-Yw[k])
        eje[k-1]=exp(dv/ae)*eje[k]*Ye[k]+(ej[k]-eji[k])*(1-Ye[k])-dv*P[k]
        eji[k-1]=exp(dv/ai)*eji[k]*Yi[k]+(ej[k]-eje[k])*(1-Yi[k])
        
        ij[k-1]=ij[k]*Yw[k]+(ije[k]+iji[k])*(1-Yw[k])
        ije[k-1]=exp(dv/ae)*ije[k]*Ye[k]+(ij[k]-iji[k])*(1-Ye[k])
        iji[k-1]=exp(dv/ai)*iji[k]*Yi[k]+(ij[k]-ije[k])*(1-Yi[k])-dv*P[k]
    
    end

    ######################################################################
    # match v<0 and v>0 solns at v=0 using excitatory and inhibitory flux
    # this fixes re&se and ri&si
    ######################################################################
    
    re=(sje[n0m]*(eji[n0p]-eji[n0m])-sji[n0m]*(eje[n0p]-eje[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m]);
    se=(rje[n0p]*(eji[n0p]-eji[n0m])-rji[n0p]*(eje[n0p]-eje[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m]);
    
    ri=(sje[n0m]*(iji[n0p]-iji[n0m])-sji[n0m]*(ije[n0p]-ije[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m]);
    si=(rje[n0p]*(iji[n0p]-iji[n0m])-rji[n0p]*(ije[n0p]-ije[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m]);

    return r,re,ri
    
end
    
#########################################################################
#########################################################################
# ThIn for modulated shot EIF with E and I 
#########################################################################
#########################################################################
function EIFShotCurrentThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,fkHz)
    
    # find the stable and unstable fixed points
    vs,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)
    
    # get useful steady-state parameters
    r,v,P,~,~=EIFShotCurrentThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT)
    
    v=collect(vlb-dv/2:dv:vth-dv/2)
    n=length(v)
    w=2*pi*fkHz
    F=(dT*exp.((v .-vT)/dT)-v)/tau

    # the reset
    kre=findall(v.>vre)[1]
    kreA=zeros(n)
    kreA[kre]=1
    
    # these are required later
    dFdvu=(exp((vu-vT)/dT)-1)/tau
    dFdvs=(exp((vs-vT)/dT)-1)/tau

    # points just above and below upper fixed point
    nup=findall(v.>vu)[1]
    num=nup-1
    
    # points just above and below lower fixed point
    nsp=findall(v.>vs)[1]
    nsm=nsp-1
    
    X=[dv*(1.0./F[1:n-1]+1.0./F[2:n])/2;0]
    
    Xw=exp.(-1im*w*X); Yw=1.0./Xw
    Xe=exp.(-Re*X);  Ye=1.0./Xe
    Xi=exp.(-Ri*X);  Yi=1.0./Xi
    
    dup=v[nup]-vu;
    Pu=P[nup]-dup*(P[nup+1]-P[nup])/dv
    
    #########################################################
    # NOW THE MODULATION
    #########################################################
   
    # intitialise the arrays for the various solutions
    (sj,sp,sje,sji)=(zeros(Complex{Float64},n) for i=1:4)
    (rj,rp,rje,rji)=(zeros(Complex{Float64},n) for i=1:4)
    (ej,ep,eje,eji)=(zeros(Complex{Float64},n) for i=1:4)
    (ij,ip,ije,iji)=(zeros(Complex{Float64},n) for i=1:4)
    
    # these are required for the >vs cases
    (pA,pB,jA,jB,jeA,jeB,jiA,jiB)=(zeros(Complex{Float64},n) for i=1:8)  # homogeneous solutions   
    (pE,pI,jE,jI,jeE,jeI,jiE,jiI)=(zeros(Complex{Float64},n) for i=1:8)  # inhomogeneous cases Re and Ri

    #########################################################
    #
    # integrate from vlb to vs 
    # vlb -----> vs ...vre... vu ..... vth
    #
    #########################################################
    
    sj[1]=0; sje[1]=0; sji[1]=-1; # the boundary condition at vlb
    ej[1]=0; eje[1]=0; eji[1]=0; 
    ij[1]=0; ije[1]=0; iji[1]=0; 
    
    for k=1:1:nsm-1
        
        sj[k+1]=sj[k]*Xw[k]+(sje[k]+sji[k])*(1-Xw[k])
        sje[k+1]=exp(-dv/ae)*(sje[k]*Xe[k]+(sj[k]-sji[k])*(1-Xe[k]))
        sji[k+1]=exp(-dv/ai)*(sji[k]*Xi[k]+(sj[k]-sje[k])*(1-Xi[k]))
        
        ej[k+1]=ej[k]*Xw[k]+(eje[k]+eji[k])*(1-Xw[k])
        eje[k+1]=exp(-dv/ae)*(eje[k]*Xe[k]+(ej[k]-eji[k])*(1-Xe[k]))+P[k]*dv
        eji[k+1]=exp(-dv/ai)*(eji[k]*Xi[k]+(ej[k]-eje[k])*(1-Xi[k]))
        
        ij[k+1]=ij[k]*Xw[k]+(ije[k]+iji[k])*(1-Xw[k])
        ije[k+1]=exp(-dv/ae)*(ije[k]*Xe[k]+(ij[k]-iji[k])*(1-Xe[k]))
        iji[k+1]=exp(-dv/ai)*(iji[k]*Xi[k]+(ij[k]-ije[k])*(1-Xi[k]))+P[k]*dv 
    
    end

    #########################################################
    #
    # now integrate from vu to vth
    # vlb .....vs ...vre... vu -----> vth
    #
    # we need to integrate two homogeneous solutions A and B where 
    # because at v_u je+ji=j we use the two initial conditions
    # jA=1 jeA=1 and jiA=0 and 
    # jB=0 jeB=1 and jiB=-1 
    #
    #########################################################
    jA[nup]= 1; jeA[nup]=1; jiA[nup]=0;
    jB[nup]= 0; jeB[nup]=1; jiB[nup]=-1
    jE[nup]= 0; jeE[nup]=0; jiE[nup]=0
    jI[nup]= 0; jeI[nup]=0; jiI[nup]=0
    
    #########################################################
    # do the run from vu - > vth
    #########################################################
    
    for k=nup:n-1
        
        jA[k+1]=jA[k]*Xw[k]+(jeA[k]+jiA[k])*(1-Xw[k])
        jeA[k+1]=exp(-dv/ae)*(jeA[k]*Xe[k]+(jA[k]-jiA[k])*(1-Xe[k]))
        jiA[k+1]=exp(-dv/ai)*(jiA[k]*Xi[k]+(jA[k]-jeA[k])*(1-Xi[k]))
            
        jB[k+1]=jB[k]*Xw[k]+(jeB[k]+jiB[k])*(1-Xw[k])
        jeB[k+1]=exp(-dv/ae)*(jeB[k]*Xe[k]+(jB[k]-jiB[k])*(1-Xe[k]))
        jiB[k+1]=exp(-dv/ai)*(jiB[k]*Xi[k]+(jB[k]-jeB[k])*(1-Xi[k]))
        
        jE[k+1]=jE[k]*Xw[k]+(jeE[k]+jiE[k])*(1-Xw[k])
        jeE[k+1]=exp(-dv/ae)*(jeE[k]*Xe[k]+(jE[k]-jiE[k])*(1-Xe[k]))+dv*P[k]
        jiE[k+1]=exp(-dv/ai)*(jiE[k]*Xi[k]+(jE[k]-jeE[k])*(1-Xi[k]))
        
        jI[k+1]=jI[k]*Xw[k]+(jeI[k]+jiI[k])*(1-Xw[k])
        jeI[k+1]=exp(-dv/ae)*(jeI[k]*Xe[k]+(jI[k]-jiI[k])*(1-Xe[k]))
        jiI[k+1]=exp(-dv/ai)*(jiI[k]*Xi[k]+(jI[k]-jeI[k])*(1-Xi[k]))+dv*P[k] 
    
    end
    
    # sort out the boundaries for the Re part
    
    a=(jiE[n]*jB[n]-jE[n]*jiB[n])/(jA[n]*jiB[n]-jiA[n]*jB[n]);
    b=(jiE[n]*jA[n]-jE[n]*jiA[n])/(jB[n]*jiA[n]-jiB[n]*jA[n]);
    
    ej[nup:n]= a*jA[nup:n] +  b*jB[nup:n] + jE[nup:n]
    ep[nup:n]= a*pA[nup:n] +  b*pB[nup:n] + pE[nup:n]
    eje[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]+ jeE[nup:n]
    eji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]+ jiE[nup:n]
    
    # sort out the boundaries for the Ri part
    
    a=(jiI[n]*jB[n]-jI[n]*jiB[n])/(jA[n]*jiB[n]-jiA[n]*jB[n]);
    b=(jiI[n]*jA[n]-jI[n]*jiA[n])/(jB[n]*jiA[n]-jiB[n]*jA[n])
    
    ij[nup:n]= a*jA[nup:n] +  b*jB[nup:n] + jI[nup:n]
    ip[nup:n]= a*pA[nup:n] +  b*pB[nup:n] + pI[nup:n]
    ije[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]+ jeI[nup:n]
    iji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]+ jiI[nup:n]

    # for the r part we require that the values at threshold are 
    # jr=1 and ji=0 so that
    
    a=jiB[n]/(jA[n]*jiB[n]-jiA[n]*jB[n])
    b=jiA[n]/(jB[n]*jiA[n]-jiB[n]*jA[n])
    
    rj[nup:n]= a*jA[nup:n] +  b*jB[nup:n]
    rp[nup:n]= a*pA[nup:n] +  b*pB[nup:n]
    rje[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]
    rji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]
    
    #########################################################
    #
    # now integrate from vu to vs
    # vlb .....vs <---vre--- vu ..... vth
    #
    #########################################################
    
    # Initial conditions
    for q in (rj,rp,rje,rji, ej,ep,eje,eji, ij,ip,ije,iji) q[num]=q[nup] end
   
    for k=num:-1:nsp+1

        rj[k-1]=rj[k]*Yw[k]+(rje[k]+rji[k])*(1-Yw[k]) - kreA[k]
        rje[k-1]=exp(dv/ae)*rje[k]*Ye[k]+(rj[k]-rji[k])*(1-Ye[k])
        rji[k-1]=exp(dv/ai)*rji[k]*Yi[k]+(rj[k]-rje[k])*(1-Yi[k])
        
        ej[k-1]=ej[k]*Yw[k]+(eje[k]+eji[k])*(1-Yw[k])
        eje[k-1]=exp(dv/ae)*eje[k]*Ye[k]+(ej[k]-eji[k])*(1-Ye[k])-dv*P[k]
        eji[k-1]=exp(dv/ai)*eji[k]*Yi[k]+(ej[k]-eje[k])*(1-Yi[k])
        
        ij[k-1]=ij[k]*Yw[k]+(ije[k]+iji[k])*(1-Yw[k])
        ije[k-1]=exp(dv/ae)*ije[k]*Ye[k]+(ij[k]-iji[k])*(1-Ye[k])
        iji[k-1]=exp(dv/ai)*iji[k]*Yi[k]+(ij[k]-ije[k])*(1-Yi[k])-dv*P[k]
    
    end

    rp[nsp]=(rj[nsp]-rje[nsp]-rji[nsp])/F[nsp]
    ep[nsp]=(ej[nsp]-eje[nsp]-eji[nsp])/F[nsp]
    ip[nsp]=(ij[nsp]-ije[nsp]-iji[nsp])/F[nsp]
    
    re=(sje[nsm]*(eji[nsp]-eji[nsm])-sji[nsm]*(eje[nsp]-eje[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);
    se=(rje[nsp]*(eji[nsp]-eji[nsm])-rji[nsp]*(eje[nsp]-eje[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);
    
    ri=(sje[nsm]*(iji[nsp]-iji[nsm])-sji[nsm]*(ije[nsp]-ije[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);
    si=(rje[nsp]*(iji[nsp]-iji[nsm])-rji[nsp]*(ije[nsp]-ije[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);

    return r,re,ri

end

#########################################################################
#########################################################################
# ThIn for modulated shot LIF with E and I 
#########################################################################
#########################################################################
function LIFShotConductThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,fkHz,ee,ei)

    # get useful steady-state parameters
    r,v,P,Je,Ji=LIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,ee,ei)
 
    v=collect(vlb+dv/2:dv:vth-dv/2)
    n=length(v)
    F=-v/tau
    
    # find points above and below the origin
    n0p=findall(v.>0)[1]
    n0m=n0p-1
    
    # the reset
    kre=minimum(findall(v.>vre))
    kreA=zeros(n)
    kreA[kre]=1
    
    # convert PSPs at v=0 to conductance parameters
    be=ae/ee
    bi=ai/ei
    betae=1/be-1
    betai=1/bi-1

    #########################################################################
    #########################################################################
    # NOW THE MODULATION
    #########################################################################
    #########################################################################
    w=2*pi*fkHz;
    X=[dv*(1.0./F[1:n-1] .+1.0./F[2:n])/2;0]

    Xw=exp.(-1im*w*X); Yw=1.0 ./Xw;
    Xe=exp.(-Re*X);    Ye=1.0 ./Xe;
    Xi=exp.(-Ri*X);    Yi=1.0 ./Xi;
    
    # intitialise the arrays for the various solutions
    (sj,sp,sje,sji)=(zeros(Complex{Float64},n) for i=1:4)
    (rj,rp,rje,rji)=(zeros(Complex{Float64},n) for i=1:4)
    (ej,ep,eje,eji)=(zeros(Complex{Float64},n) for i=1:4)
    (ij,ip,ije,iji)=(zeros(Complex{Float64},n) for i=1:4)
 
    ##############################################################
    # integrate over the range      vlb >-----> 0 ....vre.... vth
    ##############################################################
    sj[1]=0; sje[1]=0; sji[1]=-1; # the boundary condition at vlb
    ej[1]=0; eje[1]=0; eji[1]=0;  
    ij[1]=0; ije[1]=0; iji[1]=0
    
    for k=1:1:n0m-1
        
        Xbe=((ee-v[k+1])/(ee-v[k]))^betae
        Xbi=((ei-v[k+1])/(ei-v[k]))^betai
        
        sj[k+1]=sj[k]*Xw[k]+(sje[k]+sji[k])*(1-Xw[k])
        sje[k+1]=Xbe*(sje[k]*Xe[k]+(sj[k]-sji[k])*(1-Xe[k]))
        sji[k+1]=Xbi*(sji[k]*Xi[k]+(sj[k]-sje[k])*(1-Xi[k]))
        
        ej[k+1]=ej[k]*Xw[k]+(eje[k]+eji[k])*(1-Xw[k]);
        eje[k+1]=Xbe*(eje[k]*Xe[k]+(ej[k]-eji[k])*(1-Xe[k]))+P[k]*dv
        eji[k+1]=Xbi*(eji[k]*Xi[k]+(ej[k]-eje[k])*(1-Xi[k]))
        
        ij[k+1]=ij[k]*Xw[k]+(ije[k]+iji[k])*(1-Xw[k])
        ije[k+1]=Xbe*(ije[k]*Xe[k]+(ij[k]-iji[k])*(1-Xe[k]))
        iji[k+1]=Xbi*(iji[k]*Xi[k]+(ij[k]-ije[k])*(1-Xi[k]))+P[k]*dv; 

    end

    ###############################################################
    # integrate over the range      vlb ........ 0 <---vre---< vth
    ###############################################################
    rj[n]=1; rje[n]=1; rji[n]=0 # the boundary condition at vth
    ej[n]=0; eje[n]=0; eji[n]=0 
    ij[n]=0; ije[n]=0; iji[n]=0 

    for k=n:-1:n0p+1
        
        Ybe=((ee-v[k-1])/(ee-v[k]))^betae;
        Ybi=((ei-v[k-1])/(ei-v[k]))^betai
        
        rj[k-1]=rj[k]*Yw[k]+(rje[k]+rji[k])*(1-Yw[k]) - kreA[k]
        rje[k-1]=Ybe*rje[k]*Ye[k]+(rj[k]-rji[k])*(1-Ye[k])
        rji[k-1]=Ybi*rji[k]*Yi[k]+(rj[k]-rje[k])*(1-Yi[k])
        
        ej[k-1]=ej[k]*Yw[k]+(eje[k]+eji[k])*(1-Yw[k])
        eje[k-1]=Ybe*eje[k]*Ye[k]+(ej[k]-eji[k])*(1-Ye[k])-dv*P[k];
        eji[k-1]=Ybi*eji[k]*Yi[k]+(ej[k]-eje[k])*(1-Yi[k])
        
        ij[k-1]=ij[k]*Yw[k]+(ije[k]+iji[k])*(1-Yw[k])
        ije[k-1]=Ybe*ije[k]*Ye[k]+(ij[k]-iji[k])*(1-Ye[k])
        iji[k-1]=Ybi*iji[k]*Yi[k]+(ij[k]-ije[k])*(1-Yi[k])-dv*P[k]
    
    end

    ######################################################################
    # match v<0 and v>0 solns at v=0 using excitatory and inhibitory flux
    # this fixes re&se and ri&si
    ######################################################################
    
    re=(sje[n0m]*(eji[n0p]-eji[n0m])-sji[n0m]*(eje[n0p]-eje[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m])
    se=(rje[n0p]*(eji[n0p]-eji[n0m])-rji[n0p]*(eje[n0p]-eje[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m])
    
    ri=(sje[n0m]*(iji[n0p]-iji[n0m])-sji[n0m]*(ije[n0p]-ije[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m])
    si=(rje[n0p]*(iji[n0p]-iji[n0m])-rji[n0p]*(ije[n0p]-ije[n0m]))/(rje[n0p]*sji[n0m]-rji[n0p]*sje[n0m])
    
    return r,re,ri
    
end

#########################################################################
#########################################################################
# ThIn for modulated shot EIF with E and I 
#########################################################################
#########################################################################
function EIFShotConductThInEI1(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,fkHz,ee,ei)

    # find the stable and unstable fixed points
    vs,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)

    # get useful steady-state parameters
    r,v,P,Je,Ji=EIFShotConductThInEI0(dv,vth,vre,vlb,Re,ae,Ri,ai,tau,vT,dT,ee,ei)
    
    v=collect(vlb-dv/2:dv:vth-dv/2)
    n=length(v);
    w=2*pi*fkHz;
    F=(dT*exp.((v .-vT)/dT) .-v)/tau;
    
   # the reset
    kre=findall(v.>vre)[1]
    kreA=zeros(n)
    kreA[kre]=1

    # convert PSPs at v=0 to conductance parameters
    be=ae/ee;
    bi=ai/ei;
    betae=1/be-1;
    betai=1/bi-1;
    
    # these are required later
    dFdvu=(exp((vu-vT)/dT)-1)/tau;
    dFdvs=(exp((vs-vT)/dT)-1)/tau;
    
    # points just above and below upper fixed point
    nup=findall(v.>vu)[1]
    num=nup-1
    
    # points just above and below lower fixed point
    nsp=findall(v.>vs)[1]
    nsm=nsp-1

    X=[dv*(1.0./F[1:n-1]+1.0./F[2:n])/2;0]
    
    Xw=exp.(-1im*w*X); Yw=1.0./Xw
    Xe=exp.(-Re*X);  Ye=1.0./Xe
    Xi=exp.(-Ri*X);  Yi=1.0./Xi
    
    dup=v[nup]-vu;
    Pu=P[nup]-dup*(P[nup+1]-P[nup])/dv

    #########################################################
    # NOW THE MODULATION
    #########################################################

    # intitialise the arrays for the various solutions
    (sj,sp,sje,sji)=(zeros(Complex{Float64},n) for i=1:4)
    (rj,rp,rje,rji)=(zeros(Complex{Float64},n) for i=1:4)
    (ej,ep,eje,eji)=(zeros(Complex{Float64},n) for i=1:4)
    (ij,ip,ije,iji)=(zeros(Complex{Float64},n) for i=1:4)
    
    # these are required for the >vs cases
    (pA,pB,jA,jB,jeA,jeB,jiA,jiB)=(zeros(Complex{Float64},n) for i=1:8)  # homogeneous solutions   
    (pE,pI,jE,jI,jeE,jeI,jiE,jiI)=(zeros(Complex{Float64},n) for i=1:8)  # inhomogeneous cases Re and Ri 
    
    #########################################################################
    #
    # integrate from vlb to vs 
    # vlb -----> vs ...vre... vu ..... vth
    #
    #########################################################################
    sj[1]=0; sje[1]=0; sji[1]=-1; # the boundary condition at vlb
    ej[1]=0; eje[1]=0; eji[1]=0; 
    ij[1]=0; ije[1]=0; iji[1]=0; 
    
    for k=1:1:nsm-1

        Xbe=((ee-v[k+1])/(ee-v[k]))^betae
        Xbi=((ei-v[k+1])/(ei-v[k]))^betai
        
        sj[k+1]=sj[k]*Xw[k]+(sje[k]+sji[k])*(1-Xw[k])
        sje[k+1]=Xbe*(sje[k]*Xe[k]+(sj[k]-sji[k])*(1-Xe[k]))
        sji[k+1]=Xbi*(sji[k]*Xi[k]+(sj[k]-sje[k])*(1-Xi[k]))
        
        ej[k+1]=ej[k]*Xw[k]+(eje[k]+eji[k])*(1-Xw[k])
        eje[k+1]=Xbe*(eje[k]*Xe[k]+(ej[k]-eji[k])*(1-Xe[k]))+P[k]*dv
        eji[k+1]=Xbi*(eji[k]*Xi[k]+(ej[k]-eje[k])*(1-Xi[k]))
        
        ij[k+1]=ij[k]*Xw[k]+(ije[k]+iji[k])*(1-Xw[k])
        ije[k+1]=Xbe*(ije[k]*Xe[k]+(ij[k]-iji[k])*(1-Xe[k]))
        iji[k+1]=Xbi*(iji[k]*Xi[k]+(ij[k]-ije[k])*(1-Xi[k]))+P[k]*dv
    
    end
    
    #########################################################################
    #
    # now integrate from vu to vth
    # vlb .....vs ...vre... vu -----> vth
    #
    # we need to integrate two homogeneous solutions A and B where 
    # because at v_u je+ji=j we use the two initial conditions
    # jA=1 jeA=1 and jiA=0 and 
    # jB=0 jeB=1 and jiB=-1 
    #
    #########################################################################
    jA[nup]= 1; jeA[nup]=1; jiA[nup]=0;
    jB[nup]= 0; jeB[nup]=1; jiB[nup]=-1
    jE[nup]= 0; jeE[nup]=0; jiE[nup]=0
    jI[nup]= 0; jeI[nup]=0; jiI[nup]=0
    
    #########################################################
    # do the run from vu - > vth
    #########################################################
    for k=nup:n-1
    
        Xbe=((ee-v[k+1])/(ee-v[k]))^betae
        Xbi=((ei-v[k+1])/(ei-v[k]))^betai;

        jA[k+1]=jA[k]*Xw[k]+(jeA[k]+jiA[k])*(1-Xw[k]);
        jeA[k+1]=Xbe*(jeA[k]*Xe[k]+(jA[k]-jiA[k])*(1-Xe[k]));
        jiA[k+1]=Xbi*(jiA[k]*Xi[k]+(jA[k]-jeA[k])*(1-Xi[k])); 

        jB[k+1]=jB[k]*Xw[k]+(jeB[k]+jiB[k])*(1-Xw[k]);
        jeB[k+1]=Xbe*(jeB[k]*Xe[k]+(jB[k]-jiB[k])*(1-Xe[k]));
        jiB[k+1]=Xbi*(jiB[k]*Xi[k]+(jB[k]-jeB[k])*(1-Xi[k])); 
    
        jE[k+1]=jE[k]*Xw[k]+(jeE[k]+jiE[k])*(1-Xw[k]);
        jeE[k+1]=Xbe*(jeE[k]*Xe[k]+(jE[k]-jiE[k])*(1-Xe[k]))+dv*P[k];
        jiE[k+1]=Xbi*(jiE[k]*Xi[k]+(jE[k]-jeE[k])*(1-Xi[k])); 
    
        jI[k+1]=jI[k]*Xw[k]+(jeI[k]+jiI[k])*(1-Xw[k]);
        jeI[k+1]=Xbe*(jeI[k]*Xe[k]+(jI[k]-jiI[k])*(1-Xe[k]));
        jiI[k+1]=Xbi*(jiI[k]*Xi[k]+(jI[k]-jeI[k])*(1-Xi[k]))+dv*P[k]; 
   
    end

    # sort out the boundaries for the Re part
    
    a=(jiE[n]*jB[n]-jE[n]*jiB[n])/(jA[n]*jiB[n]-jiA[n]*jB[n]);
    b=(jiE[n]*jA[n]-jE[n]*jiA[n])/(jB[n]*jiA[n]-jiB[n]*jA[n]);
    
    ej[nup:n]= a*jA[nup:n] +  b*jB[nup:n] + jE[nup:n]
    ep[nup:n]= a*pA[nup:n] +  b*pB[nup:n] + pE[nup:n]
    eje[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]+ jeE[nup:n]
    eji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]+ jiE[nup:n]
    
    # sort out the boundaries for the Ri part
    
    a=(jiI[n]*jB[n]-jI[n]*jiB[n])/(jA[n]*jiB[n]-jiA[n]*jB[n]);
    b=(jiI[n]*jA[n]-jI[n]*jiA[n])/(jB[n]*jiA[n]-jiB[n]*jA[n])
    
    ij[nup:n]= a*jA[nup:n] +  b*jB[nup:n] + jI[nup:n]
    ip[nup:n]= a*pA[nup:n] +  b*pB[nup:n] + pI[nup:n]
    ije[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]+ jeI[nup:n]
    iji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]+ jiI[nup:n]

    # for the r part we require that the values at threshold are 
    # jr=1 and ji=0 so that
    
    a=jiB[n]/(jA[n]*jiB[n]-jiA[n]*jB[n])
    b=jiA[n]/(jB[n]*jiA[n]-jiB[n]*jA[n])
    
    rj[nup:n]= a*jA[nup:n] +  b*jB[nup:n]
    rp[nup:n]= a*pA[nup:n] +  b*pB[nup:n]
    rje[nup:n]=a*jeA[nup:n] + b*jeB[nup:n]
    rji[nup:n]=a*jiA[nup:n] + b*jiB[nup:n]

    #########################################################################
    #
    # now integrate from vu to vs
    # vlb .....vs <---vre--- vu ..... vth
    #
    #########################################################################

    # Initial conditions
    for q in (rj,rp,rje,rji, ej,ep,eje,eji, ij,ip,ije,iji) q[num]=q[nup] end

    for k=num:-1:nsp+1

        Ybe=((ee-v[k-1])/(ee-v[k]))^betae;
        Ybi=((ei-v[k-1])/(ei-v[k]))^betai;

        rj[k-1]=rj[k]*Yw[k]+(rje[k]+rji[k])*(1-Yw[k]) - kreA[k]
        rje[k-1]=Ybe*rje[k]*Ye[k]+(rj[k]-rji[k])*(1-Ye[k]);
        rji[k-1]=Ybi*rji[k]*Yi[k]+(rj[k]-rje[k])*(1-Yi[k]);
    
        ej[k-1]=ej[k]*Yw[k]+(eje[k]+eji[k])*(1-Yw[k]);
        eje[k-1]=Ybe*eje[k]*Ye[k]+(ej[k]-eji[k])*(1-Ye[k])-dv*P[k];
        eji[k-1]=Ybi*eji[k]*Yi[k]+(ej[k]-eje[k])*(1-Yi[k]);

        ij[k-1]=ij[k]*Yw[k]+(ije[k]+iji[k])*(1-Yw[k]);
        ije[k-1]=Ybe*ije[k]*Ye[k]+(ij[k]-iji[k])*(1-Ye[k]);
        iji[k-1]=Ybi*iji[k]*Yi[k]+(ij[k]-ije[k])*(1-Yi[k])-dv*P[k];

    end
    
    rp[nsp]=(rj[nsp]-rje[nsp]-rji[nsp])/F[nsp];
    ep[nsp]=(ej[nsp]-eje[nsp]-eji[nsp])/F[nsp];
    ip[nsp]=(ij[nsp]-ije[nsp]-iji[nsp])/F[nsp];

    re=(sje[nsm]*(eji[nsp]-eji[nsm])-sji[nsm]*(eje[nsp]-eje[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);
    se=(rje[nsp]*(eji[nsp]-eji[nsm])-rji[nsp]*(eje[nsp]-eje[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);

    ri=(sje[nsm]*(iji[nsp]-iji[nsm])-sji[nsm]*(ije[nsp]-ije[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);
    si=(rje[nsp]*(iji[nsp]-iji[nsm])-rji[nsp]*(ije[nsp]-ije[nsm]))/(rje[nsp]*sji[nsm]-rji[nsp]*sje[nsm]);

    return r,re,ri
 
end

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# SIMULATION FUNCTIONS
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# LIF simulation in response to time-dependent stimulation by Ret & Rit
#########################################################################
#########################################################################
function LIFShotCurrentSimEI(tau,vth,vre,vlb,t,Ret,Rit,ae,ai,v1)

    dt=t[2]-t[1]
    nt=length(t)
    
    v=ones(nt); v[1]=v1
    r=zeros(nt)
    
    # Allow for constant-rate arguments
    length(Ret)==1 ? Ret=Ret*ones(nt) : ()
    length(Rit)==1 ? Rit=Rit*ones(nt) : ()
    
    se=(rand(nt).<Ret*dt).*rand(Exponential(ae),nt) 
    si=-(rand(nt).<Rit*dt).*rand(Exponential(abs(ai)),nt) 
    
    for k=1:nt-1
        v[k+1]=v[k]-(dt/tau)*v[k]+se[k]+si[k]
        v[k+1]>vth ? (v[k+1]=vre; v[k]=vth; r[k]=1/dt) : ()
        v[k+1]<vlb ?  v[k+1]=vlb : ()
    end
    
    vend=v[end]
    rsum=sum(r)/nt
    
    return v,vend,r,rsum

end

#########################################################################
#########################################################################
# LIF simulation in response to time-dependent stimulation by Ret & Rit
#########################################################################
#########################################################################
function LIFShotConductSimEI(tau,vth,vre,vlb,t,Ret,Rit,ae,ai,ee,ei,v1)
    
    # for the simulation
    dt=t[2]-t[1]
    nt=length(t)

    # conductance parameters
    be=ae/ee
    bi=ai/ei
    betae=1/be-1;
    betai=1/bi-1;
    
    v=ones(nt); v[1]=v1
    r=zeros(nt)
    
    # Allow for constant-rate arguments
    length(Ret)==1 ? Ret=Ret*ones(nt) : ()
    length(Rit)==1 ? Rit=Rit*ones(nt) : ()
    
    se=(rand(nt).<Ret*dt).*rand(Beta(1,betae),nt) 
    si=(rand(nt).<Rit*dt).*rand(Beta(1,betai),nt) 
    
    # deals with LIF case
    dT==0 ? vth=vT : ()
    
    for k=1:nt-1
        v[k+1]=v[k]-(dt/tau)*v[k]+(ee-v[k])*se[k]+(ei-v[k])*si[k]
        v[k+1]>vth ? (v[k+1]=vre; v[k]=vth; r[k]=1/dt) : ()
        v[k+1]<vlb ?  v[k+1]=vlb : ()
    end
    
    vend=v[end]
    rsum=sum(r)/nt
    
    return v,vend,r,rsum

end



#########################################################################
#########################################################################
# Simulation in response to time-dependent stimulation by Ret & Rit
# NB With dT=0 this simulates an LIF with vth=vT
#########################################################################
#########################################################################
function EIFShotCurrentSimEI(tau,dT,vT,vth,vre,vlb,t,Ret,Rit,ae,ai,v1)

    dt=t[2]-t[1]
    nt=length(t)
    
    v=ones(nt); v[1]=v1
    r=zeros(nt)
    
    # Allow for constant-rate arguments
    length(Ret)==1 ? Ret=Ret*ones(nt) : ()
    length(Rit)==1 ? Rit=Rit*ones(nt) : ()
    
    se=(rand(nt).<Ret*dt).*rand(Exponential(ae),nt) 
    si=-(rand(nt).<Rit*dt).*rand(Exponential(abs(ai)),nt) 
    
    # deals with LIF case. 
    dT==0 ? vth=vT : ()
    
    for k=1:nt-1
        v[k+1]=v[k]+(dt/tau)*(dT*exp((v[k]-vT)/dT)-v[k])+se[k]+si[k]
        v[k+1]>vth ? (v[k+1]=vre; v[k]=vth; r[k]=1/dt) : ()
        v[k+1]<vlb ?  v[k+1]=vlb : ()
    end
    
    vend=v[end]
    rsum=sum(r)/nt
    
    return v,vend,r,rsum

end

#########################################################################
#########################################################################
# Simulation in response to time-dependent stimulation by Ret & Rit
# NB With dT=0 this simulates an LIF with vth=vT
#########################################################################
#########################################################################
function EIFShotConductSimEI(tau,dT,vT,vth,vre,vlb,t,Ret,Rit,ae,ai,ee,ei,v1)
    
    # for the simulation
    dt=t[2]-t[1]
    nt=length(t)

    # conductance parameters
    be=ae/ee
    bi=ai/ei
    betae=1/be-1;
    betai=1/bi-1;
    
    v=ones(nt); v[1]=v1
    r=zeros(nt)
    
    # Allow for constant-rate arguments
    length(Ret)==1 ? Ret=Ret*ones(nt) : ()
    length(Rit)==1 ? Rit=Rit*ones(nt) : ()
    
    se=(rand(nt).<Ret*dt).*rand(Beta(1,betae),nt) 
    si=(rand(nt).<Rit*dt).*rand(Beta(1,betai),nt) 
    
    # deals with LIF case
    dT==0 ? vth=vT : ()
    
    for k=1:nt-1
        v[k+1]=v[k]+(dt/tau)*(dT*exp((v[k]-vT)/dT)-v[k])+(ee-v[k])*se[k]+(ei-v[k])*si[k]
        v[k+1]>vth ? (v[k+1]=vre; v[k]=vth; r[k]=1/dt) : ()
        v[k+1]<vlb ?  v[k+1]=vlb : ()
    end
    
    vend=v[end]
    rsum=sum(r)/nt
    
    return v,vend,r,rsum

end
