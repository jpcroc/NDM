#!/usr/bin/python3
# -*-coding:Utf-8 -*
import json
import math
import sys

import numpy as np
import scipy as sp
from scipy import optimize
import matplotlib.pyplot as plt

KB = 8.6173303E-5  # constante de boltzmann en eV/K

#X0_EFERMI = 1.0
EFLOOP = np.arange(0., 2.9, 0.1)



class Phase:
	def __init__(self, **defaut_attributes):    
		for attr_name, attr_value in defaut_attributes.items():
			setattr(self, attr_name, attr_value)    
		self.site=[]
		self.defect=[]
		self.NA0=()
		self.NA=[]
		self.mu=[]
		self.muS=[]
		self.Pv=[]
		self.Nc=[]
		self.NAtot0=0
		print ("name",self.name)
		with open(f"{self.name}.out", "w") as ph_out :
			print("open ", self.name)
		
	def print_cdef(self):
		for defect in self.defect:
			print("defprt",defect.name)
			defect.print_c()
		return
	def write_cdef(self,var,temp):
		print ("nom",self.name)
		stch=self.stoich()
		invstch=1/stch
		asd=abs(invstch-2)
		muaw=self.mu[0]
		mubw=self.mu[1]
		with open(f"{self.name}.{temp}.out", "a") as ph_out:
			ph_out.write(f"   {var} {invstch} {asd} {muaw}  {mubw} ")
		for defect in self.defect:
			print("defwrt",defect.name)
			with open(f"c_{defect.name}.{temp}.out", "a") as def_out:
				def_out.write(f"   {var} ")
			defect.write_c()
			with open(f"{self.name}.{temp}.out", "a") as ph_out:
				ph_out.write(f"   {defect.c_tot} ")
		with open(f"{self.name}.{temp}.out", "a") as ph_out:
			ph_out.write(f"\n")

		return
		
	def set_NA(self,lprt=False):
		self.NA=self.NA0[:]
		da=[]
		for iat in range(len(self.NA)):
			da.append(0.)

		# print("preccalc",self.NA)
		for defect in self.defect:
			isit=defect.site
			nbs=self.site[isit].nb
			# print("defect site nbs ",defect.name, isit,nbs)		
			

			for iat in range(len(self.NA)):
				da[iat]=da[iat]+defect.dAB[iat]*defect.c_tot*nbs
		for iat in range(len(self.NA)):
			self.NA[iat]=self.NA[iat]+da[iat]
			stch=self.NA[1]/self.NA[0]
		if lprt :
			print("set_NA",self.species,self.NA,stch)
	def stoich(self):
		self.set_NA()
		stp=self.NA[1]/self.NA[0]
		return(stp)
	
	def calcNcPv(self,temp):
		beta=1 / (KB*temp)
		Nc=0.
		Pv=0.
		nelec=0.
		for Ee in self.dos :
			if Ee <= self.EBV :
				Pv=Pv+np.exp((Ee-self.EBV)*beta)
				nelec=nelec+1.0
			elif Ee >= self.EBC :
				Nc=Nc+np.exp((self.EBC-Ee)*beta)
			else :
				print("WTF NCPV")
				exit()
		# print(nelec)
		self.Nc.append(Nc)
		self.Pv.append(Pv)

	
class Defect:
	def __init__(self, **defaut_attributes):
		for attr_name, attr_value in defaut_attributes.items():
			setattr(self, attr_name, attr_value)
		# c is a list of concentrations per charge
		self.c = []  #q: 0 for q in self.q}
		# c_tot total concentration et1of defect with all charges (H_i+ + H_i-)
		self.c_tot = 0
		# with open(f"c_{self.name}.{temp}.out", "w") as def_out :
			# print("open ", self.name)

	def set_c(self, efermi,mu):
		# set concentration depending on mu and efermi
		# print("BETA",BETA)
		Et1=0
		Et2=0
		self.c_tot=0
		# print("mu",mu)
		# print("dab",self.dAB)
		# for atyp in range(len(self.dAB)):
			# Et1=Et1+self.dAB[atyp]*mu[atyp]
		for dab,mua in zip(self.dAB,mu):
			Et1=Et1-dab*mua		
		for q in self.q:
			iq=self.q.index(q)
			Et=Et1+self.Er[iq]+q*efermi
			# # print("CALC",self.name,self.q[iq],Et,Et1,q*efermi,self.Er[iq])
			# # print("MULT",self.name,iq,self.q[iq],self.mult[iq],self.prefactor[iq])
			# # print("EXPOS",self.name,q,Et)
			try:
				self.c[iq] = self.prefactor[iq] *self.mult[iq] * math.exp(
					-(Et * BETA))
			except OverflowError:
				self.c[iq] = self.prefactor[iq] * self.mult[iq]*np.exp(
					-(Et * BETA))
			# print("CONC",self.name,self.q[iq],self.c[iq])
					
			self.c_tot=self.c_tot+self.c[iq]
		return
	def print_c(self):
		for q in self.q:
			iq=self.q.index(q)
			print("CONC",self.name,self.q[iq],self.c[iq])
		print("CONCTOT",self.name,self.c_tot)
		print(" ")
	def write_c(self):
		with open(f"c_{self.name}.{temp}.out", "a") as def_out:
			def_out.write(f" {self.c_tot}")
			for q in self.q:
				iq=self.q.index(q)
				def_out.write(f"   {self.c[iq]}")
			def_out.write(f"\n")



class Site:
	def __init__(self, **defaut_attributes):
		# Attributes name, charge, energy and reaction to the defect
		for attr_name, attr_value in defaut_attributes.items():
			setattr(self, attr_name, attr_value)

def phase_creation(INP,Ph):
	PHR=[]
	"""Constructs list of objects Default and phaseeters."""
	with open(INP, 'r') as inp:
		for line in inp:
			next_line = inp.readline()
			# Creation of PHASE 
			if line.split()[0] == "Phase":
				while len(next_line.strip()) != 0:
					if next_line[0] != '#':
						PHR.append(
							Phase(**json.loads(next_line.rstrip("\n\r ")))
						)
					next_line = inp.readline()
			# Creation of DEFECTS list
			elif line.split()[0] == "Sites":
				while len(next_line.strip()) != 0:
					if next_line[0] != '#':
						PHR[0].site.append(
							Site(**json.loads(next_line.rstrip("\n\r ")))
						)
					next_line = inp.readline()
			elif line.split()[0] == "Defects":
				while len(next_line.strip()) != 0:
					print(next_line)
					if next_line[0] != '#':
						PHR[0].defect.append(
							Defect(**json.loads(next_line.rstrip("\n\r ")))
						)
					next_line = inp.readline()
		# PHR[0].NA0=[]
		PHR[0].efermi=0.5*(PHR[0].EBC-PHR[0].EBV)
		for speci in range(len(PHR[0].species)):
			PHR[0].NA.append(0.)
		for siti in range(len(PHR[0].site)):
			try:
				typ=PHR[0].site[siti].typ
				PHR[0].NA[typ]=PHR[0].NA[typ]+PHR[0].site[siti].nb
				print("Site",PHR[0].site[siti].name,"nb ",PHR[0].site[siti].nb,"type ",PHR[0].site[siti].typ)
			except:
				print("Site",PHR[0].site[siti].name,"nb ",PHR[0].site[siti].nb,"empty site")
		# for speci in range(len(PHR[0].species)):
			# #print("na0",PHR[0].NA0[spec])
			# print("sp",PHR[0].species[speci],PHR[0].NA0[speci])
		# PHR[0].NA0=(PHR[0].NA0[:])
		print("PERFECT STOICH",PHR[0].species,PHR[0].NA)
		PHR[0].NA0=(PHR[0].NA)
		PHR[0].NAtot0=sum(PHR[0].NA[:])
		Eperat=PHR[0].E0SC/(PHR[0].NSC*PHR[0].NAtot0)
		PHR[0].E0=PHR[0].E0SC/PHR[0].NSC
		# for speci in range(len(PHR[0].species)):
			# PHR[0].mu.append(Eperat)
			# PHR[0].muS.append(Eperat)
			# print("bounds",PHR[0].species[speci],PHR[0].bounds[speci])
		PHR[0].mu.append(float(PHR[0].muAR))
		PHR[0].muS.append(float(PHR[0].muAR))
		EperU=(PHR[0].E0+float(PHR[0].muAR)*PHR[0].NA[0])/PHR[0].NA[1]
		PHR[0].mu.append(EperU)
		PHR[0].muS.append(EperU)

		
		print(" ")
		print("DEFECTS")
		for defect in PHR[0].defect:
			try:
				print(defect.name,"site",defect.site,"delta AB", defect.dAB)
			except:
				print(defect.name, "no dAB really ?")
				exit()
			try:
				print("prefactor",defect.name,defect.prefactor)#
			except:
				defect.prefactor=[]
				for q in defect.q:
					defect.prefactor.append(1.0)
				print("prefactor",defect.name,defect.prefactor, "by default")#doStuff(a.property)
			try:
				print("FVibT",defect.name,defect.FVibT)#doStuff(a.property)
			except:
				defect.FVibT=[]
				for temp in PHR[0].TEMP:
					defect.FVibT.append(0.0)
				print("FVibT",defect.name,defect.FVibT, "by default")#doStuff(a.property)
			try:
				print("mult",defect.name,defect.mult)#doStuff(a.property)
			except:
				defect.mult=[]
				for q in defect.q:
					defect.mult.append(1)
				print("mult",defect.name,defect.mult, "by default")#doStuff(a.property)
			defect.Er=[]
			for q in defect.q:
				iq=defect.q.index(q)
				Er=defect.energy[iq]-PHR[0].E0SC+q*PHR[0].EBV
				defect.Er.append(Er)
				defect.c.append(0.0)
			for q,er in zip(defect.q,defect.Er):
				print("Er",defect.name,q,er)#doStuff(a.property)
		
			
			
			print(" ")
	# print(PHR[0].dosf)		
	if hasattr(PHR[0],"dosf") :
		PHR[0].dos = np.genfromtxt(PHR[0].dosf)
		for temp in PHR[0].TEMP:
			print(temp)
			PHR[0].calcNcPv(temp)
		print(PHR[0].name,".Nc",PHR[0].Nc)
		print(PHR[0].name,".Pv",PHR[0].Pv)
	Ph.append(PHR[0])	
	del PHR
			
			
def stoichT(x,Ph):
	Ph.set_NA()
	st=Ph.NA0[0]*Ph.NA[1]-((Ph.NA0[1]+x)*Ph.NA[0])
	return(st)

def stoichzero(muA,x,Ph):
	Ph.mu[0]=muA
	Ph.mu[1]=(Ph.E0-muA*Ph.NA0[0])/Ph.NA0[1]
	print("muA_muB",Ph.mu[0],Ph.mu[1])
	ef0=Ph.efermi
	ef1 = float(sp.optimize.fsolve(charge_totale, ef0, (Ph)))
	stz=stoichT(x,Ph)
	return(stz)

def F_mu2(muA,muB,Ph) :
	frE=FreeE(Ph)
	Ph.set_NA()
	fmu=frE-muA*Ph.NA[0]-muB*Ph.NA[1]
	fr0=muA*Ph.NA[0]+muB*Ph.NA[1]
	# print("fmu ",fmu,"fr0 ",fr0,"frE", frE)
	fmu2=fmu**2
	return(fmu2)

def stoich2 (x,Ph):
	st2=(stoichT(x,Ph))**2
	return (st2)

		
def FreeE(Ph,lprt=False) :
	eint=Eint(Ph)
	mts=mTS(Ph)
	F=Eint(Ph)+mTS(Ph)
	if lprt:
		print("FreeE", F, "EINT ", eint, "-TS ",mts)
	
	return(F)

def Eint(Ph):
	eint=Ph.E0
	# print("EINT base", eint)
	for defect in Ph.defect :
		isit=defect.site
		nbs=Ph.site[isit].nb
		for q in defect.q :
			iq=defect.q.index(q)
			dEdq=nbs*defect.c[iq]*defect.Er[iq]
			eint=eint+dEdq
	# print("EINT def", eint)
	return(eint)

def decalFvib(Ph,itemp):
	for defect in Ph.defect:
		# print (itemp,"def Erav",defect.name, defect.Er,"fvib",defect.FVibT[itemp])
		for ier in range(len(defect.Er)) :
			er=defect.Er[ier]+defect.FVibT[itemp]
			defect.Er[ier]=er
		# print ("def Erap",defect.name, defect.Er)
		# defect.Er=list(map(lambda x:x+defect.FVibT[itemp], defect.Er.copy()))
		# print ("def ErXX",defect.name, defect.Er)			

def DedecalFvib(Ph,itemp):
	for defect in Ph.defect:
		defect.Er=list(map(lambda x:x-defect.FVibT[itemp], defect.Er.copy()))
	

def mTS(Ph):
	mts=0.
	CS0=[]
	dcs0=[]
	for site in Ph.site:
		CS0.append(1.0)
		dcs0.append(0.0)
	for defect in Ph.defect :
		isit=defect.site
		nbs=Ph.site[isit].nb
		for q in defect.q :
			iq=defect.q.index(q)
			dmts=nbs*defect.c[iq]*math.log(defect.c[iq]/defect.mult[iq])/BETA
			mts=mts+dmts
			dcs0[isit]=dcs0[isit]+defect.c[iq]
	for site in Ph.site:
		isite=Ph.site.index(site)
		if hasattr(site,'typ'):
			CS0[isite]=CS0[isite]-dcs0[isite]
			# print("CS0 dcs0",isite,CS0[isite],dcs0[isite])
			dmts=site.nb*CS0[isite]*math.log(CS0[isite])/BETA
			mts=mts+dmts
		else:
			CS0[isite]=0.
	# print("mtS ",mts)
	return(mts)
	
	
def charge_totale(efermi,Ph):
	"Calculate the total charge of material, with the equation \sum_{q, D}q[D_q]."""
	q_tot=0
	for defect in Ph.defect:
		defect.set_c(efermi,Ph.mu)
		isit=defect.site
		nbs=Ph.site[isit].nb
		# print("def",defect.name,defect.site,nbs,defect.q,defect.c)
		for ql,cl in zip(defect.q,defect.c):
			# print("Q C",ql,cl)
			q_tot=q_tot+ql * cl*nbs
	# print("BETA",efermi,BETA,q_tot)

	if hasattr(Ph,"dosf"):
		itemp=Ph.TEMP.index(temp)
		print("BETA",efermi,BETA,itemp,Ph.Pv[itemp],Ph.Nc[itemp])

		q_tot=q_tot+Ph.Pv[itemp]*np.exp(BETA*(-efermi))-Ph.Nc[itemp]*np.exp(BETA*(Ph.EBV+efermi-Ph.EBC))
		Ph.nh=Ph.Pv[itemp]*np.exp(BETA*(-efermi))
		Ph.ne=Ph.Nc[itemp]*np.exp(BETA*(Ph.EBV+efermi-Ph.EBC))
		print("ELH",Ph.ne,Ph.nh)
	# q_tot = sum([ql * cl for defect in Ph.defect for ql,cl in zip(defect.c,defect.q)])
	# Calculation of electrons and holes concentrations
	# q_tot += PV * math.exp(-efermi * BETA)
	# q_tot -= NC * math.exp(-(EBC - (efermi+EBV)) * BETA)
	# print("q_tot = {q_tot}")
	# print ("cht",efermi,q_tot)
	return(q_tot)

def	Fmintot_mu(mu,*args):
	# global xt
	# print(" ")
	# print("MUTEST",mu)
	# x=xt
	Ph=args[0]
	x=args[1]
	Ph.mu[0]=mu[0]
	Ph.mu[1]=mu[1]
	ef0=Ph.efermi
	ef1=float(sp.optimize.fsolve(charge_totale, ef0, (Ph)))
	fmu2=F_mu2(Ph.mu[0],Ph.mu[1],Ph)
	sta=abs(stoichT(x,Ph))
	fmt=fmu2+sta
	# try:
		# rati=fmu2/sta
		# # print("FMT",fmt,"FMU2",fmu2,"sta",sta,"ratio",rati)
	# except ZeroDivisionError:
		# # print ("You can't divide by zero!")
		
	return(fmt)	

def	Fmintot_x(xmuB,*args):
	# global xt
	# print(" ")
	# print("MUTEST",mu)
	# x=xt
	Ph=args[0]
	
	Ph.mu[1]=xmuB[1]
	x=xmuB[0]
	ef0=Ph.efermi
	ef1=float(sp.optimize.fsolve(charge_totale, ef0, (Ph)))
	fmu2=F_mu2(Ph.mu[0],Ph.mu[1],Ph)
	sta=abs(stoichT(x,Ph))
	fmt=fmu2+sta
	# try:
		# rati=fmu2/sta
		# # print("FMT",fmt,"FMU2",fmu2,"sta",sta,"ratio",rati)
	# except ZeroDivisionError:
		# # print ("You can't divide by zero!")
		
	return(fmt)	

def calcmu(Ph) :
	# global PHT
	print(" ")
	print("**********CALCMU**********")
	print(" ")

	Ph.set_NA()
	print("******SIMPLE SOLUTION*********")

	muAt=Ph.mu[0]
	print("muAT",muAt)

	muA=float(sp.optimize.fsolve(stoichzero, muAt, (xt,Ph)))
	muB=(Ph.E0-muA*Ph.NA0[0])/Ph.NA0[1]
	Ph.muS[0]=Ph.mu[0]
	Ph.muS[1]=Ph.mu[1]

	Ph.mu[0]=muA
	Ph.mu[1]=muB
	# Etot=Ph.mu[0]*Ph.NA0[0]+Ph.mu[1]*Ph.NA0[1]
	efermi = float(sp.optimize.fsolve(charge_totale, X0_EFERMI, (Ph)))
	Ph.efermi=efermi
	print("EFermi charge_totale", Ph.efermi,charge_totale(Ph.efermi,Ph))
	if hasattr(Ph,"dosf") :
		print("nelec",Ph.ne)
		print("nholes",Ph.nh)

	Ph.set_NA(lprt=True)
	Ph.print_cdef()

	Ftot=FreeE(Ph,lprt=True) #Eint(Ph)+mTS(Ph)
	print("FTOTmu0",Ftot,Ph.E0,Ph.mu[0],Ph.mu[1])

	print(" ")
	print("******REAL SOLUTION*********")

	# exit()
	# bounds=[Ph.bounds[ispec] for ispec in range(len(Ph.species))]
	# print(bounds)
	delta=0.5
	muam=Ph.mu[0]-delta
	muaM=Ph.mu[0]+delta
	mubm=Ph.mu[1]-delta
	mubM=Ph.mu[1]+delta

	bounds=[(muam,muaM),(mubm,mubM)]
	print("bounds",bounds)
	critere, securite = 1, 0
	while critere > 1E-10 and securite < 20:
		print("Essai numéro ", securite)
		result = sp.optimize.differential_evolution(
			Fmintot_mu, bounds=bounds, args=(Ph,xt),
			strategy='best2bin', maxiter=1000, popsize=10, tol=1E-15,
			mutation=(0.1, 0.5), recombination=0.7, disp=True,#updating='deferred'
		)
		critere = result.fun
		securite += 1
		print(critere, securite,result.success,result.x)
		
	Ph.mu[0]=result.x[0]
	Ph.mu[1]=result.x[1]
	efermi = float(sp.optimize.fsolve(charge_totale, X0_EFERMI, (Ph)))
	Ph.efermi=efermi
	print("EFermi charge_totale", Ph.efermi,charge_totale(Ph.efermi,Ph))
	Ph.set_NA(lprt=True)
	Ph.print_cdef()
	Ftot=FreeE(Ph,lprt=True) #Eint(Ph)+mTS(Ph)
	print("FTOTmu0",Ftot,Ph.E0,Ph.mu[0],Ph.mu[1])
	print("Dmu",Ph.mu[0]-Ph.muS[0],Ph.mu[1]-Ph.muS[1])


def calcx(Ph,simple=False) :
	global muA
	print(" ")
	print("*******CALCX************")
	print(" ")
	# print("mu0",Ph.mu[0],Ph.mu[1])

	Ph.set_NA()
	# print("mu",Ph.mu[0],Ph.mu[1])

	# print(" ")
	print("******SIMPLE SOLUTION*********")

	Ph.mu[1]=(Ph.E0-Ph.mu[0]*Ph.NA0[0])/Ph.NA0[1]
	print("mu",Ph.mu[0],Ph.mu[1])

	efermi = float(sp.optimize.fsolve(charge_totale, X0_EFERMI, (Ph)))
	Ph.efermi=efermi
	print("EFermi charge_totale", Ph.efermi,charge_totale(Ph.efermi,Ph))
	Ph.set_NA(lprt=True)
	Ph.print_cdef()

	Ftot=FreeE(Ph,lprt=True) #Eint(Ph)+mTS(Ph)
	print("CHECK",muA,Ph.mu[0])
	st0=Ph.stoich()
	print("FTOTmu0",Ftot,Ph.E0,Ph.mu[0],Ph.mu[1])

	if simple :
		Ph.set_NA(lprt=True)
		Ph.print_cdef()
		Ftot=FreeE(Ph,lprt=True) #Eint(Ph)+mTS(Ph)
		return
	
	Ph.muS[0]=Ph.mu[0]
	Ph.muS[1]=Ph.mu[1]

	print(" ")
	print("******REAL SOLUTION*********")
	delta=0.01
	mubm=Ph.mu[1]-delta
	mubM=Ph.mu[1]+delta

	bounds=[(-0.05,0.1),(mubm,mubM)]
	print("bounds calcX",bounds)
	critere, securite = 1, 0
	while critere > 1E-10 and securite < 20:
		print("Essai numéro ", securite)
		result = sp.optimize.differential_evolution(
			Fmintot_x, bounds=bounds, args=(Ph,),
			strategy='best2bin', maxiter=4000, popsize=10, tol=1E-15,
			mutation=(0.1, 0.5), recombination=0.7, disp=True,#updating='deferred'
		)
		critere = result.fun
		securite += 1
		print(critere, securite,result.success,result.x)
	print("CHECK",muA,Ph.mu[0])
		
	#Ph.mu[0]=result.x[0]
	Ph.mu[1]=result.x[1]
	efermi = float(sp.optimize.fsolve(charge_totale, X0_EFERMI, (Ph)))
	Ph.efermi=efermi
	print("EFermi charge_totale", Ph.efermi,charge_totale(Ph.efermi,Ph))
	Ph.set_NA(lprt=True)
	Ph.print_cdef()
	Ftot=FreeE(Ph,lprt=True) #Eint(Ph)+mTS(Ph)
	print("FTOTmu0",Ftot,Ph.E0,Ph.mu[0],Ph.mu[1])
	st1=Ph.stoich()
	print("Dmu2",Ph.mu[0],Ph.mu[0]-Ph.muS[0],Ph.mu[1]-Ph.muS[1],st0,st1)


PHT=[]
phase_creation("inpT1",PHT)
X0_EFERMI = PHT[0].efermi


for itemp in range(len(PHT[0].TEMP)):
	temp=PHT[0].TEMP[itemp]
	BETA=1 / (KB*temp)

	xt=0.0
	decalFvib(PHT[0],itemp)
	calcmu(PHT[0])
	if hasattr(PHT[0],"dosf") :
		print("nelec",PHT[0].ne)
		print("nholes",PHT[0].nh)

	print("mu1",PHT[0].mu[0],PHT[0].mu[1])
	with open(f"muA.out", "w") as muout :
		print("open MU")
	muA0=PHT[0].mu[0]
	muB0=PHT[0].mu[1]
	mumin=PHT[0].mu[0]-0.
	mumax=PHT[0].mu[0]+3.5
	
	MULOOP=np.arange(mumin,mumax, 0.02 )
	for muA in MULOOP:
		with open(f"muA.{temp}.out", "a") as muout :
			muout.write(f"   {muA} ")
			
		PHT[0].mu[0]=muA
		# dmuA=muA-muA0
		# PHT[0].mu[1]=muB0-dmuA*PHT[0].NA0[0]/PHT[0].NA0[1]
	
		calcx(PHT[0])
		# calcx(PHT[0],simple=True)
	
		PHT[0].write_cdef(muA,temp)
		stch=PHT[0].stoich()
		invstch1=1/stch
		Ftot=FreeE(PHT[0])
		with open(f"muA.{temp}.out", "a") as muout :
			muout.write(f"   {Ftot} {PHT[0].mu[0]} {PHT[0].mu[1]} {stch} {invstch1} \n")
	DedecalFvib(PHT[0],itemp)

