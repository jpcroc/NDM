#!/usr/bin/python3
import numpy as np
import sys
class page :
	def __init__(self,ip) :
		self.ip=ip
		self.text=[]
		self.sals=[]
class Sal: 
	def __init__(self,isal,itranche,indpage) :
		self.diplomes=[]
		self.isal=isal
		self.tranche=itranche
		self.page=indpage
		self.comp=True
		self.lok=True
		self.av=[]
		self.affect=''
		for il in range(0,10,1) :		
			self.av.append(" ")
			
	def check(self,testA,remp):
		self.lok=True
		if not hasattr(self,testA):
			self.lok=False
			print('SAL sans ', testA, 'isal tranche', self.isal,self.tranche)
			self.comp=False
			if testA=='nom' :
				tranchepb=self.tranche
				for textc in Pages[ipg].text :
					if textc.col==0 and textc.tranche==tranchepb: 
						if textc.x==35. and 'INGEN' in textc.text:
							self.nom=textc.text.replace('INGENIEUR','')
							self.comp=True
							self.lok=True
				if not self.lok	:
					exit()
			if testA=='affect' :
				self.affect=remp
				self.comp=True
				self.lok=True
			if testA=='ech' :
				self.affect=remp
				self.comp=True
				self.lok=True
			if testA=='embauche' :
				self.embauche='pas trouvé'
				self.comp=True
				self.lok=True
			if testA=='clembauche' :
				self.clembauche='pas trouvé'
				self.comp=True
				self.lok=True
		if testA=='ech' :
			if int(self.ech)==1 or int(self.ech)==2 or int(self.ech)==3 :
				self.cat=1
			if int(self.ech)==4 or int(self.ech)==5 :
				self.cat=2
			if int(self.ech)==6 or int(self.ech)==7  :
				self.cat=3
		if not self.lok :
			print('problem page tranche',testA, self.page, self.tranche)
			print(self.__dict__)
			exit()


class textel:
	def __init__(self,x,y,xf,yf,ipg,text,itxt):
		self.text=text.replace("'","")
		self.x=x
		if 'ADMINISTRATIF' in self.text :
			self.text=text.replace("CADRE ADMINISTRATIF","")
			self.x=146.

		self.y=y
		self.xf=xf
		self.yf=yf
		self.page=ipg
		self.itxt=itxt
	def ftranche(self) :
		self.tranche=-1
		y=self.y
		yf=self.yf
		trlmin=[]
		trlmax=[]
		trlmin.append(80.)
		trlmax.append(141.)
		trlmin.append(135.)
		trlmax.append(197.)
		trlmin.append(191.)
		trlmax.append(253.)
		trlmin.append(247.)
		trlmax.append(310.)
		trlmin.append(303.)
		trlmax.append(365.)
		trlmin.append(360.)
		trlmax.append(418.)
		for il in range(0,6,1) :
			if yf<trlmax[il] and y> trlmin[il] :
				self.tranche=il
				# print("text tranche", self.tranche,self.y,self.yf,self.text)
				break
		if self.tranche==-1 :
			if y> 78 and yf< 420 :
				print("text PAS DE tranche", self.y,self.yf,self.text,self.x,self.xf)
				exit()
		
	def fcol(self) :
		self.col=-1		
		x=self.x
		xf=self.xf
		trcmin=[]
		trcmax=[]
		trcmin.append(30.)
		trcmax.append(147.)
		trcmin.append(145.)
		trcmax.append(182.)
		trcmin.append(182.)
		trcmax.append(374.)
		trcmin.append(374.)
		trcmax.append(528.)
		trcmin.append(528.)
		trcmax.append(811.)
		if self.tranche >-1 :
			for il in range(0,5,1) :
				if x<trcmax[il] and x> trcmin[il] :
					self.col=il
					break
			if self.col==-1:
				if x> 30.0 and xf< 811. :
					print("text PAS DE col", self.x,self.xf)
					exit()
		
nf=sys.argv[1]
with open(nf) as file:
# file = open('outT2.N2')
	lines=file.readlines()
# file.close()
#print(lines)
nfo=nf+".csv"
fileo=open(nfo,'w')
ipg=-1

itextdec=0
itxt=0
Pages=[]
Textlist=[]
tranchelue=[]
tranche2sal=[]
Sallist=[]
for itr in range (0,6,1):
	tranchelue.append('0')
	tranche2sal.append('-1')
ldpg=[]

for iline in range(len(lines)):
	line=lines[iline]
	# print(iline,line)
	elem=line.split()
	# print(line)
	# print(elem )
	if 'Line' in elem[0] :
		itextdec=0
		# print('LINE')
	else :
		itextdec=itextdec+1
		if itextdec==1 :
			ipg=ipg+1
			# print('nouvpage ',ipg,iline)
			ldpg.append(iline)
	lignefin=iline

print(*ldpg)

for ipg	in range(len(ldpg)):
	ipd=ldpg[ipg]
	try:
		ipf=ldpg[ipg+1]-1
	except :
		ipf=lignefin
	print(ipg,ipd,ipf)
	print("nouvpage", ipg,ipd,ipf)
	Pages.append(page(ipg))	
	pagecur=Pages[ipg]
	pasfini=True
	# while pasfini :
		# écrire la ligne dans un fichier
		
	# do while pas bon
		# lire les lignes du fichier
		
	nsal=-1
	for itr in range (0,6,1):
		tranchelue[itr]=0
		tranche2sal[itr]=-1

	with open(nf) as file:
		lines=file.readlines()[ipd:ipf]
		for iline in range(len(lines)):
			indline=iline+ipd
			line=lines[iline]
			# print(iline,line)
			elem=line.split()
			# print(line)
			# print(elem )
			if 'Line' in elem[0] :
				itextdec=0
				# print('LINE')
			else :
				# itextdec=itextdec+1
				# if itextdec==1 :
					# ipg=ipg+1
					# print('nouvpage ',ipg,iline)
					# Pages.append(page(ipg))
					# itxt=-1		
				# print(" ")
				# print('TEXT')
				elemextr=elem.copy()
				# print ('elem',elem)
				del elemextr[0:5]
				# print ('extr',elemextr)
				# for el in elemextr :
					# elemextr=elemextr+""
				mrg=" ".join(elemextr)
				# print ("mrg",mrg)
				x=float(elem[1])
				y=float(elem[2])
				xf=float(elem[3])
				yf=float(elem[4])
				textT=textel(x,y,xf,yf,ipg,mrg,itxt)
				textT.ftranche()
				# print(textT.text,textT.x)
				itranche=textT.tranche
				if itranche >-1:
					if tranchelue[itranche]==0 :
						tranchelue[itranche]=1
						nsal=nsal+1
						pagecur.sals.append(Sal(nsal,itranche,ipg))
						tranche2sal[textT.tranche]=nsal
					salcur=pagecur.sals[tranche2sal[itranche]]
					textT.fcol()
					if textT.col>-1 :
						itxt=itxt+1
						Pages[ipg].text.append(textT)
						# textc=Pages[ipg].text[itxt]
					# textc.ftranche()
					# textc.fcol()
						Textlist.append(textT)
						
					if textT.col==0 :
						if textT.x>58 and textT.x<60:
							salcur.nom=textT.text
						if textT.x>113 and textT.x<115 :
							salcur.age=textT.text
						if textT.x>33 and textT.x<35 and textT.xf>73 and textT.xf< 75 :
							salcur.badge=textT.text
						if textT.x>33 and textT.x<35 and textT.xf< 45 :
							salcur.genre='H'
						if textT.x>33 and textT.x<35 and textT.xf>45 and textT.xf< 60 :
							salcur.genre='F'

					if textT.col==1 :
						if textT.y>360 and textT.y<366 :
							salcur.affect=salcur.affect+textT.text
						if textT.y>304 and textT.y<310 :
							salcur.affect=salcur.affect+textT.text
						if textT.y>248 and textT.y<254 :
							salcur.affect=salcur.affect+textT.text
						if textT.y>192 and textT.y<197 :
							salcur.affect=salcur.affect+textT.text
						if textT.y>136 and textT.y<141 :
							salcur.affect=salcur.affect+textT.text
						if textT.y>80 and textT.y<85 :
							salcur.affect=salcur.affect+textT.text
					if textT.col==2 :
						salcur.diplomes.append(textT.text)
					if textT.col==3 :
						if textT.x>375 and textT.x <377:
							if textT.y>122 and textT.y<127 :
								salcur.embauche=textT.text
							if textT.y>178 and textT.y<183 :
								salcur.embauche=textT.text
							if textT.y>234 and textT.y<239 :
								salcur.embauche=textT.text
							if textT.y>290 and textT.y<295 :
								salcur.embauche=textT.text
							if textT.y>346 and textT.y<351 :
								salcur.embauche=textT.text
							if textT.y>402 and textT.y<407 :
								salcur.embauche=textT.text
						if textT.x>448 and textT.x <450:
							if textT.y>122 and textT.y<127 :
								salcur.clembauche=textT.text
							if textT.y>178 and textT.y<183 :
								salcur.clembauche=textT.text
							if textT.y>234 and textT.y<239 :
								salcur.clembauche=textT.text
							if textT.y>290 and textT.y<295 :
								salcur.clembauche=textT.text
							if textT.y>346 and textT.y<351 :
								salcur.clembauche=textT.text
							if textT.y>402 and textT.y<407 :
								salcur.clembauche=textT.text
						if textT.x>473 and textT.x <475:
							if textT.y>365 and textT.y<370 :
								salcur.classement=textT.text
							if textT.y>309 and textT.y<314 :
								salcur.classement=textT.text
							if textT.y>252 and textT.y<258 :
								salcur.classement=textT.text
							if textT.y>195 and textT.y<202 :
								salcur.classement=textT.text
							if textT.y>139 and textT.y<146 :
								salcur.classement=textT.text
							if textT.y>85 and textT.y<90 :
								salcur.classement=textT.text								
						if textT.x>375 and textT.x <377:
							if textT.y>365 and textT.y<370 :
								salcur.ech=textT.text[-2:]
							if textT.y>309 and textT.y<314 :
								salcur.ech=textT.text[-2:]
							if textT.y>252 and textT.y<258 :
								salcur.ech=textT.text[-2:]
							if textT.y>195 and textT.y<202 :
								salcur.ech=textT.text[-2:]
							if textT.y>139 and textT.y<146 :
								salcur.ech=textT.text[-2:]
							if textT.y>85 and textT.y<90 :
								salcur.ech=textT.text[-2:]								

					if textT.col==4 :
						# print("AV",textT.x,textT.text)
						unaff=0
						if textT.x>532 and textT.x <534:
							salcur.av[0]=textT.text
							unaff=1
						if textT.x>560 and textT.x <565:
							salcur.av[1]=textT.text
							unaff=1
						if textT.x>588 and textT.x <590:
							salcur.av[2]=textT.text
							unaff=1
						if textT.x>616 and textT.x <618:
							salcur.av[3]=textT.text
							unaff=1
						if textT.x>644 and textT.x <646:
							salcur.av[4]=textT.text
							unaff=1
						if textT.x>671 and textT.x <673:
							salcur.av[5]=textT.text
							unaff=1
						if textT.x>699 and textT.x <701:
							salcur.av[6]=textT.text
							unaff=1
						if textT.x>727 and textT.x <729:
							salcur.av[7]=textT.text
							unaff=1
						if textT.x>755 and textT.x <758:
							salcur.av[8]=textT.text
							unaff=1
						if textT.x>783 and textT.x <785:
							salcur.av[9]=textT.text
							unaff=1
						if unaff==0:
							print("UNAFF AVANCEMENT",textT.x,textT.text)
							exit()
						
						
		for isc in range(len(pagecur.sals)):
			salcur=pagecur.sals[isc]
			salcur.diplome=" ".join(salcur.diplomes)
			delattr(salcur,'diplomes')
		# for salcur in pagecur.sals:
			salcur.check('nom','')
			salcur.check('age','')
			salcur.check('badge','')
			salcur.check('ech','')
			remp=''
			try:
				remp=pagecur.sals[isc-1].affect
			except :
				pass
			try:
				remp=pagecur.sals[isc-1].affect
			except :
				pass
			# except :
				# remp=pagecur.sals[isc+1].affect
			salcur.check('affect',remp)
			salcur.check('embauche','')
			salcur.check('clembauche','')			 
		# print("pass2")
		for salcur in pagecur.sals:
			
			# salcur.check('affect')
			if not salcur.comp :
				print('PROBLEM')
				exit()
			# print(" ")
			# print(salcur.__dict__)
			# print(salcur.badge,salcur.genre)
			# print(" ")
			# print(salcur.ech)
			fileo.write(f"{salcur.badge} ; {salcur.nom}; {salcur.genre}; {salcur.age}; {salcur.cat}; {salcur.ech}; {salcur.classement}; {salcur.affect}; {salcur.diplome}; {salcur.embauche};{salcur.clembauche} ; {salcur.av[0]}; {salcur.av[1]}; {salcur.av[2]}; {salcur.av[3]}; {salcur.av[4]}; {salcur.av[5]}; {salcur.av[6]}; {salcur.av[7]}; {salcur.av[8]}; {salcur.av[9]} \n")
			# for cl in salcur.av:
				# fileo.write(f"{cl} ;")
				# fileo.write(f" \n")
# for itr in range (0,6,1):
	# for ipg	in range(len(ldpg)):
# for ipg	in range(len(ldpg)):
	# for itr in range (0,6,1):
		# for text in Textlist :
			# # print (text)
			# if text.page==ipg and text.col ==3 and text.tranche==itr:
				# print (text.page, text.tranche, text.x,text.y,text.xf,text.yf,text.text)
