#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

CONSTANT NumContings=1
CONSTANT NumArms=6,SeqLenMax=6
StrCONSTANT ContingList="1;"
CONSTANT ResampleNumber=10000

function runCBAS(anInfoWv)
	wave anInfoWv

	getAllInfo()
	
	make /o/n=(1,3) critLengthInfo
	critLengthInfo[0][0]=0
	critLengthInfo[0][1]=400
	critLengthInfo[0][2]=800
	doAllSeqAcrossAn(critLengthInfo)
	
	wave wv=$"wtsHipReorder"
	string list=nameOfwave(wv)+";"
	doResampleAndFindK(list,1)
	saveExperiment
end

function getAllInfo()
	wave /t expInfo
	string list="sex;types;lesion;implant;"
	string DF
	variable i,j
	for(i=0;i<itemsInList(list);i+=1)
		make /o/n=0 outWv
		if(i==0)
			make /o/n=0 experiment
		endif
		for(j=0;j<numpnts(expInfo);j+=1)
			DF="root:"+expInfo[j]+":"
			wave thisOne=$DF+stringFromList(i,list)
			duplicate /o thisOne hold
			concatenate /NP=0 "hold;",outWv
			if(i==0)
				hold=j
				concatenate /NP=0 "hold;",experiment
			endif
		endfor
		duplicate /o outWv $stringFromList(i,list)
	endfor
	killwaves /z hold,outWv
	getContributions()
	getAllContingLength()
end

function getAnNum(which)
	variable which
	wave experiment
	variable val=which
	variable anNum=0
	if(which==0)
		return anNum
	endif
	if(which>=numpnts(experiment))
		return -1
	endif
	do
		if(experiment[val-1]==experiment[val])
			anNum+=1
		else
			return anNum
		endif
		val-=1
	while(val>0)
	return anNum
end

function getContributions()
	wave experiment
	make /o/n=(numpnts(experiment),NumContings+1) allAnContrib=0
	wave /t expInfo
	variable val
	variable i
	for(i=0;i<numpnts(experiment);i+=1)
		val=getAnNum(i)
		wave wv=$"root:"+expInfo[experiment[i]]+":reAllOneAn"+num2str(val)
		getContingInfo(wv)
		wave contingInfo
		allAnContrib[i][0,min(dimsize(contingInfo,0)-1,NumContings)]=1
		if(i==0)
			duplicate /o/r=[][0,1] contingInfo allContingInfo
		else
			if(combineContingInfo()==-1)
				print i
				print "Contingencies not aligned"
			endif
		endif
	endfor
end

function getContingInfo(wv)
	wave wv
	getCotingNum(wv)
	wave numTrials,wConting,isExp
	make /o/n=(numpnts(wConting)) hold=1
	if(numpnts(hold)>1)
		hold[1,numpnts(hold)-1]=wConting[p]-wConting[p-1]
	endif
	hold=selectnumber(hold[p]==1,NaN,p)
	wavetransform zapNaNs,hold
	make /o/n=(numpnts(hold),3) contingInfo
	variable val
	variable i
	for(i=0;i<numpnts(hold);i+=1)
		if(i==0)
			val=0
		else
			val=sum(numTrials,0,hold[i]-1)
		endif
		if(isExp[hold[i]])
			contingInfo[i][0,1]=wv[val][q+3]+10
		else
			contingInfo[i][0,1]=wv[val][q+3]
		endif
		if(i<numpnts(hold)-1)
			contingInfo[i][2]=hold[i+1]-hold[i]
		else
			contingInfo[i][2]=numpnts(numTrials)-hold[i]
		endif
	endfor
end

function combineContingInfo()
	wave contingInfo,allContingInfo
	variable val=min(dimsize(allContingInfo,0),dimsize(contingInfo,0))
	make /o/n=(val,dimsize(allContingInfo,1)) hold
	hold=abs(allContingInfo[p][q]-contingInfo[p][q])
	wavestats /q hold
	if(v_max>0)
		return -1
	endif
	if(dimsize(allContingInfo,0)<dimsize(contingInfo,0))
		duplicate /o/r=[][0,1] contingInfo allContingInfo
	endif
end

function getAllContingLength()
	wave experiment
	wave /t expInfo
	wave allAnContrib
	make /o/n=(numpnts(experiment),NumContings+1) allContingLength=NaN,allContingRew=NaN
	variable val
	variable start,finish
	variable i,j
	for(i=0;i<numpnts(experiment);i+=1)
		val=getAnNum(i)
		wave wv=$"root:"+expInfo[experiment[i]]+":reAllOneAn"+num2str(val)
		getCotingNum(wv)
		wave numTrials,wConting
		make /o/n=(numpnts(wConting)) hold
		for(j=0;j<=NumContings;j+=1)
			if(allAnContrib[i][j]==1)
				hold=selectnumber(wConting[p]==j,NaN,p)
				wavestats /q hold
				allContingLength[i][j]=sum(numTrials,v_min,v_max)
				if(j==0)
					start=0
				else
					start=sum(numTrials,0,v_min-1)
				endif
				finish=sum(numTrials,0,v_max)-1
				wavestats /q/rmd=[start,finish][2] wv
				allContingRew[i][j]=v_sum
			endif
		endfor
	endfor
end

function getCotingNum(wv)
	wave wv
	wavestats /q/rmd=[][0] wv
	make /o/n=(dimsize(wv,0)) hold
	variable mn=v_min
	variable mx=v_max
	make /o/n=(mx-mn+1) numTrials,wConting,isExp=0
	variable start,finish
	variable homeS=NaN,outerS=NaN
	variable homeF=NaN,outerF=NaN
	variable i,j
	for(i=mn;i<=mx;i+=1)
		hold=selectnumber(wv[p][0]==i,NaN,p)
		wavestats /q hold
		numTrials[i-mn]=v_npnts
		start=v_min
		finish=v_max
		if(i==mn)
			wConting[i-mn]=0
			homeS=wv[start][3]
			outerS=wv[start][4]
			homeF=wv[finish][3]
			outerF=wv[finish][4]
			if(homeS>=0 && homeF>=0)
				isExp[i-mn]=0
			else
				isExp[i-mn]=1
			endif
//Experiment specific correction
		
//End experiment specific correction			
		elseif(wv[start][3]>=0)
				if(homeF!=wv[start][3] || outerF!=wv[start][4])
					wConting[i-mn]=wConting[i-mn-1]+1
					homeS=wv[start][3]
					outerS=wv[start][4]
					homeF=wv[finish][3]
					outerF=wv[finish][4]
					if(homeS>=0 && homeF>=0)
						isExp[i-mn]=0
					else
						isExp[i-mn]=1
					endif
				elseif(wv[finish][3]!=wv[start][3] || wv[finish][4]!=wv[start][4])
					wConting[i-mn]=wConting[i-mn-1]+1
					homeS=wv[start][3]
					outerS=wv[start][4]
					homeF=wv[finish][3]
					outerF=wv[finish][4]
					if(homeS>=0 && homeF>=0)
						isExp[i-mn]=0
					else
						isExp[i-mn]=1
					endif
				else
					wConting[i-mn]=wConting[i-mn-1]
					isExp[i-mn]=0
				endif
		elseif((wv[start][3]>=0)==0)
			wConting[i-mn]=wConting[i-mn-1]
			homeS=wv[start][3]
			outerS=wv[start][4]
			homeF=wv[finish][3]
			outerF=wv[finish][4]
			isExp[i-mn]=1
		endif
	endfor
	fixReExp()
end

function fixReExp()
	wave numTrials,wConting,isExp
	variable i
	for(i=1;i<numpnts(isExp)-1;i+=1)
		if(isExp[i-1]==1 && isExp[i]==0 && isExp[i+1]==1)
			isExp[i]=1
			wConting[i+1,numpnts(wConting)-1]-=1
		endif 
	endfor
	for(i=0;i<numpnts(isExp)-1;i+=1)
		if(isExp[i+1]==1 && isExp[i]==0 && wConting[i+1]==wConting[i])
			isExp[i]=1
			wConting[i,numpnts(wConting)-1]+=1
		endif 
	endfor
end

function doAllSeqAcrossAn(critLengthWv)
	wave critLengthWv
	variable howMany=SeqLenMax
	killDataFolder /z root:seqDF
	newDataFolder /o root:seqDF
	variable timer
	variable i,j,k
	for(j=1;j<=howMany;j+=1)
		printf "----Starting length %g----",j
		timer=startMStimer
		for(i=0;i<itemsInList(ContingList);i+=1)
			getManySeqAcrossAn(str2num(stringFromList(i,ContingList)),j,dimsize(critLengthWv,0))
		endfor
		printf "Finished length %g in %g\r",j,stopMStimer(timer)/1e6
	endfor
	saveExperiment
	for(k=0;k<dimsize(critLengthWv,0);k+=1)
		print "----getting crit trials----"
		timer=startMStimer
		if(critLengthWv[k][0]==0)
			getZerothOrder(critLengthWv)
		else
			getAllTrialToPerfect(critLengthWv[k][0],critLengthWv[k][1])
		endif
		wave allTrialToPerfect
		duplicate /o allTrialToPerfect $"allTrialToPerfect"+num2str(k)
		printf "Finished crit %g in %g\r",k,stopMStimer(timer)/1e6
	endfor
	saveExperiment
	print "\r----On to sequence counts----"
	for(k=0;k<dimsize(critLengthWv,0);k+=1)
		for(j=1;j<=howMany;j+=1)
			printf "----Starting length %g----\r",j
			for(i=0;i<itemsInList(ContingList);i+=1)
				timer=startMStimer
				wave seqTrialWv=$"root:seqDF:allSeqAllAn_"+stringFromList(i,ContingList)+"_"+num2str(j)
				wave critWv=$"allTrialToPerfect"+num2str(k)
				wave contribWv=allAnContrib
				wave seqCntsWv=$"root:seqDF:seqCnts_"+stringFromList(i,ContingList)+"_"+num2str(j)
				extractCritCnts(seqCntsWv,seqTrialWv,contribWv,critWv,str2num(stringFromList(i,ContingList)),k)
				printf "Finished conting %g in %g\r",str2num(stringFromList(i,ContingList)),stopMStimer(timer)/1e6
				saveExperiment
			endfor
		endfor
	endfor
	saveExperiment
end

function HowManySequences(howMany)
	variable howMany
	variable val
	variable i,j
	for(j=1;j<=howMany;j+=1)
		for(i=0;i<itemsInList(ContingList);i+=1)
			wave wv=$"root:seqDF:allSeq_"+stringFromList(i,ContingList)+"_"+num2str(j)
			val+=dimsize(wv,1)
		endfor
	endfor
	return val
end

function getManySeqAcrossAn(conting,howMany,numCrit)
	variable conting,howMany,numCrit
	getAllAnSeq(conting,howMany)
	wave seqsWv=allSequence
	wave experiment
	make /o/n=(numpnts(experiment),dimsize(seqsWv,1),numCrit) seqCnts=NaN
	variable i,j
	for(i=0;i<numpnts(experiment);i+=1)
		if(identAllSeq(seqsWv,i,conting)==1)
			wave anSeqTrial
			if(j==0)
				duplicate /o anSeqTrial allSeqAllAn
			else
				concatenate /NP=0 "anSeqTrial;",allSeqAllAn
			endif
			j+=1
			extractTotCnts(seqCnts,anSeqTrial,i)
		endif
	endfor
	make /o/n=(dimsize(allSeqAllAn,0)) hold0,hold1,hold2
	hold0=allSeqAllAn[p][0]
	hold1=allSeqAllAn[p][1]
	hold2=allSeqAllAn[p][2]
	SortColumns keyWaves={hold2,hold0,hold1},sortWaves=allSeqAllAn
	duplicate /o allSeqAllAn $"root:seqDF:allSeqAllAn_"+num2str(conting)+"_"+num2str(howMany)
	duplicate /o seqsWv $"root:seqDF:allSeq_"+num2str(conting)+"_"+num2str(howMany)
	duplicate /o seqCnts $"root:seqDF:seqCnts_"+num2str(conting)+"_"+num2str(howMany)
end

function identAllSeq(seqHold,which,conting)
	wave seqHold
	variable which,conting
	variable howMany=dimsize(seqHold,1)
	wave experiment
	wave /t expInfo
	variable val
	val=getAnNum(which)
	if(val<0)
		return -1
	endif
	wave wv=$"root:"+expInfo[experiment[which]]+":reAllOneAn"+num2str(val)
	getCotingNum(wv)
	wave numTrials,wConting
	make /o/n=(numpnts(wConting)) ses=selectnumber(wConting[p]==conting,NaN,p)
	wavetransform zapNaNs,ses
	if(numpnts(ses)==0)
		return 0
	endif
	variable i
	for(i=0;i<numpnts(ses);i+=1)
		make /o/n=(dimsize(wv,0)) oneSes
		oneSes=selectnumber(wv[p][0]==ses[i],NaN,wv[p][1]+wv[p][2]*NumArms)
		waveTransform zapNaNs,oneSes
		identSesSeq(oneSes,seqHold)
		wave seqNums
		if(i==0)
			duplicate /o seqNums anSeq
		else
			concatenate /NP=0 "seqNums;",anSeq
		endif
	endfor
	if(numpnts(anSeq)==0)
		return 0
	endif
	make /o/n=(numpnts(anSeq)) anTrial
	anTrial=selectnumber(anSeq[p]>=0,NaN,p)
	sort anSeq,anSeq,anTrial
	wavestats /q anTrial
	deletepoints v_npnts,v_numNaNs,seqNums,anTrial
	if(numpnts(anTrial)==0)
		return 0
	endif
	wavestats /q anTrial
	variable maxNum=v_max
	wavestats /q anSeq
	maxNum=max(v_max,maxNum)
	if(maxNum>(2^32) || which>(2^32))
		print "NB not enough bits in wave!!!!!!"
	endif
	make /o/i/u/n=(numpnts(anTrial),3) anSeqTrial
	anSeqTrial[][0]=which
	anSeqTrial[][1]=anTrial[p] 
	anSeqTrial[][2]=anSeq[p] 
	return 1
end

function identSesSeq(wv,seqHold)
	wave wv,seqHold
	variable howMany=dimsize(seqHold,0)
	make /o/n=(howMany) seq
	make /o/n=(numpnts(wv)) seqNums=NaN
	variable i
	for(i=0;i<numpnts(wv)-howMany+1;i+=1)
		seq=wv[p+i]
		seqNums[i]=getSequenceNum(seqHold,seq)
	endfor
end

function getSequenceNum(seqHold,seqCompWv)
	wave seqHold,seqCompWv
	duplicate /o seqCompWv hold
	variable i
	do
		if(i>=dimsize(seqHold,1))
			return -1
			print "Error: missing sequence"
		endif
		hold=abs(seqCompWv[p]-seqHold[p][i])
		if(sum(hold)==0)
			return i
		endif
		i+=1
	while(1)
end

function extractTotCnts(seqCntsWv,anSeqTrialWv,which)
	wave seqCntsWv,anSeqTrialWv
	variable which
	variable i
	make /o/n=(dimsize(anSeqTrialWv,0)) pnts=anSeqTrialWv[p][2]
	sort pnts,pnts
	for(i=numpnts(pnts)-1;i>0;i-=1)
		if(pnts[i]==pnts[i-1])
			deletepoints i,1,pnts
		endif
	endfor
	seqCntsWv[which][][0]=0
	for(i=0;i<numpnts(pnts);i+=1)
		make /o/n=(dimsize(anSeqTrialWv,0)) hold=selectnumber(anSeqTrialWv[p][2]==pnts[i],0,1)
		seqCntsWv[which][pnts[i]]=sum(hold)
	endfor
end

function extractCritCnts(seqCntsWv,seqTrialWv,contribWv,critWv,conting,critNum)
	wave seqCntsWv,seqTrialWv,contribWv,critWv
	variable conting,critNum
	variable i,val,which
	
	Variable numAn=dimsize(seqCntsWv,0)
	Variable numSeq=dimsize(seqCntsWv,1)
	Variable nthreads= ThreadProcessorCount
	Variable threadGroupID= ThreadGroupCreate(nthreads)
	
	for(which=0; which<numAn; which+=1)
		if(contribWv[which][conting]==0)
			seqCntsWv[which][][critNum]=NaN
		else
		
		make /o/n=(dimsize(seqTrialWv,0)) pnts=selectnumber(seqTrialWv[p][0]==which,NaN,p)
		wavetransform zapNaNs,pnts
		
		if(numpnts(pnts)==0)
			seqCntsWv[which][][critNum]=NaN
		else
			
			make /o/n=(numpnts(pnts),2) anSeqTrial=seqTrialWv[pnts[p]][q+1]
			
				for(val=0; val<numSeq;)
					for(i=0; i<nthreads; i+=1)
						ThreadStart threadGroupID,i,MyCounterFunc(seqCntsWv,anSeqTrial,critWv,conting,which,val,critNum)
						val+=1
						if( val>=numSeq )
							break
						endif
					endfor
					
					do
						Variable threadGroupStatus = ThreadGroupWait(threadGroupID,100)
					while( threadGroupStatus != 0 )
				endfor
			endif
		endif
	endfor
	Variable dummy= ThreadGroupRelease(threadGroupID)
end

ThreadSafe function MyCounterFunc(seqCntsWv,seqTrialWv,critWv,conting,whichAn,val,critNum)
	wave seqCntsWv,seqTrialWv,critWv
	variable conting,whichAn,val,critNum
	if(seqCntsWv[whichAn][val][0]>0)
		make /o/n=(dimsize(seqTrialWv,0)) hold=selectnumber(seqTrialWv[p][0]<=critWv[whichAn][conting] && seqTrialWv[p][1]==val,0,1)
		seqCntsWv[whichAn][val][critNum]=sum(hold)
	else
		seqCntsWv[whichAn][val][critNum]=0
	endif
	
	return stopMSTimer(-2)		// time when we finished
End

function getAllAnSeq(conting,howMany)
	variable conting,howMany
	wave experiment
	make /o/n=(0) allCnt
	make /o/b/u/n=(howMany,0) allSequence
	variable val
	variable i
	for(i=0;i<numpnts(experiment);i+=1)
		if(extractAllSeq(i,conting,howMany)==1)
			wave cnt,sequence
			getAllUnique(allSequence,allCnt,cnt,sequence)
		endif
	endfor
end

function getAllUnique(seqHold,cntHold,cntWv,seqWv)
	wave seqHold,cntHold,cntWv,seqWv
	make /o/n=(dimsize(seqWv,0)) hold
	variable i,j,k
	for(i=0;i<dimsize(seqWv,1);i+=1)
		j=0
		do
			if(j==(numpnts(cntHold)))
				make /o/b/u/n=(dimsize(seqHold,0),numpnts(cntHold)+1) $getWavesDataFolder(seqHold,2)
				seqHold[][dimsize(seqHold,1)-1]=seqWv[p][i]
				make /o/n=(numpnts(cntHold)+1) $getWavesDataFolder(cntHold,2)
				cntHold[numpnts(cntHold)-1]=cntWv[i]
				k=1
				break
			endif
			hold=abs(seqWv[p][i]-seqHold[p][j])
			if(sum(hold)==0)	
				cntHold[j]+=cntWv[i]
				break
			else
				j+=1
			endif
		while(1)
	endfor
	if(k)
		make /o/b/u/n=(dimsize(seqHold,1),dimsize(seqHold,0)) holdSeq=seqHold[q][p]
		SortColumns /r keyWaves=cntHold,sortWaves={holdSeq,cntHold}
		seqHold=holdSeq[q][p]
	endif
end

function extractAllSeq(which,conting,howMany)
	variable which,conting,howMany
	wave experiment
	wave /t expInfo
	variable val
	val=getAnNum(which)
	if(val<0)
		return -1
	endif
	wave wv=$"root:"+expInfo[experiment[which]]+":reAllOneAn"+num2str(val)
	getCotingNum(wv)
	wave numTrials,wConting
	make /o/n=(numpnts(wConting)) ses=selectnumber(wConting[p]==conting,NaN,p)
	wavetransform zapNaNs,ses
	if(numpnts(ses)==0)
		return 0
	endif
	make /o/n=0 cnt
	make /o/n=(howMany,0) sequence
	variable i,j
	for(i=0;i<numpnts(ses);i+=1)
		make /o/n=(dimsize(wv,0)) oneSes
		oneSes=selectnumber(wv[p][0]==ses[i],NaN,wv[p][1]+wv[p][2]*NumArms)
		waveTransform zapNaNs,oneSes
		if(numpnts(oneSes)>=howMany)
			extractSeq(oneSes,howMany)
			wave allSeq
			getUnique(sequence,cnt,allSeq)
			j=1
		endif
	endfor
	if(j==1)
		return 1
	else
		return 0
	endif
end

function extractSeq(wv,howMany)
	wave wv
	variable howMany
	make /o/n=(howMany) oneSeq
	variable i
	for(i=0;i<numpnts(wv)-howMany+1;i+=1)
		oneSeq=wv[p+i]
		if(i==0)
			make /o/n=(numpnts(oneSeq),1) allSeq=oneSeq
//			duplicate /o oneSeq allSeq
		else
			concatenate /NP=1 "oneSeq;",allSeq
		endif
	endfor
end

function getUnique(seqHold,cntWv,seqWv)
	wave seqHold,cntWv,seqWv
	make /o/n=(dimsize(seqWv,0)) hold
	variable i,j
	for(i=0;i<dimsize(seqWv,1);i+=1)
		j=0
		do
			if(j==(numpnts(cntWv)))
				make /o/n=(dimsize(seqHold,0),numpnts(cntWv)+1) $getWavesDataFolder(seqHold,2)
				seqHold[][dimsize(seqHold,1)-1]=seqWv[p][i]
				make /o/n=(numpnts(cntWv)+1) $getWavesDataFolder(cntWv,2)
				cntWv[numpnts(cntWv)-1]=1
				break
			endif
			hold=abs(seqWv[p][i]-seqHold[p][j])
			if(sum(hold)==0)	
				cntWv[j]+=1
				break
			else
				j+=1
			endif
		while(1)
	endfor
end

function getPerfectPerformanceSeq(wv)
	wave wv
	make /o/n=(dimsize(wv,0)) hold
	variable i,j
	for(i=0;i<dimsize(wv,1);i+=1)
		hold=floor(wv[p][i]/NumArms)
		if(sum(hold)==numpnts(hold))
			make /o/n=(j+1) perfectSeq
			perfectSeq[j]=i
			j+=1
		endif
	endfor
end

function getAnNumPerfect(conting,howMany)
	variable conting,howMany
	wave wv=$"root:seqDF:allSeqAllAn_"+num2str(conting)+"_"+num2str(howMany)
	duplicate /o wv allSeqAllAn
	wave wv=$"root:seqDF:allSeq_"+num2str(conting)+"_"+num2str(howMany)
	duplicate /o wv allSeq
	getPerfectPerformanceSeq(allSeq)
	wave perfectSeq
	wave allAnContrib
	make /o/n=(dimsize(allAnContrib,0)) anPerfect=NaN
	make /o/n=(dimsize(allSeqAllAn,0)) hold
	variable i,j
	for(i=0;i<numpnts(anPerfect);i+=1)
		hold=0
		if(allAnContrib[i][conting]==1)
			for(j=0;j<numpnts(perfectSeq);j+=1)
				hold=selectnumber(allSeqAllAn[p][0]==i && allSeqAllAn[p][2]==perfectSeq[j],hold[p],1)
			endfor
			anPerfect[i]=sum(hold)
		endif
	endfor
end

function getZerothOrder(wv)
	wave wv
	wave experiment
	wave allContingLength
	make /o/n=(numpnts(experiment),dimsize(wv,1)-1) allTrialToPerfect=NaN
	variable which
	variable i
	for(i=0;i<itemsInList(ContingList);i+=1)
		which=str2num(stringFromList(i,ContingList))
		make /o/n=(numpnts(experiment)) trialToPerfect=selectnumber(allContingLength[p][which]>=wv[0][which+1],inf,wv[0][which+1])
		allTrialToPerfect[][str2num(stringFromList(i,ContingList))]=trialToPerfect[p]
	endfor
end

function getAllTrialToPerfect(howMany,thr)
	variable howMany,thr
	wave experiment
	make /o/n=(numpnts(experiment),NumContings+1) allTrialToPerfect=NaN
	variable i
	for(i=0;i<itemsInList(ContingList);i+=1)
		if(str2num(stringFromList(i,ContingList))==0 || str2num(stringFromList(i,ContingList))==8 || str2num(stringFromList(i,ContingList))==9)
			make /o/n=(numpnts(experiment)) trialToPerfect=inf
		else
			getTrialToPerect(str2num(stringFromList(i,ContingList)),howMany,thr)
			wave trialToPerfect
		endif
		allTrialToPerfect[][str2num(stringFromList(i,ContingList))]=trialToPerfect[p]
	endfor
end

function getTrialToPerect(conting,howMany,thr)
	variable conting,howMany,thr
	wave allAnContrib
	wave wv=$"root:seqDF:allSeqAllAn_"+num2str(conting)+"_"+num2str(howMany)
	if(waveExists(wv))
		duplicate /o wv allSeqAllAn
		wave wv=$"root:seqDF:allSeq_"+num2str(conting)+"_"+num2str(howMany)
		duplicate /o wv allSeq
		getPerfectPerformanceSeq(allSeq)
		wave perfectSeq
		make /o/n=(dimsize(allAnContrib,0)) trialToPerfect=NaN
		variable i,j
		for(i=0;i<numpnts(trialToPerfect);i+=1)
			make /o/n=(dimsize(allSeqAllAn,0)) hold=NaN
			if(allAnContrib[i][conting]==1)
				for(j=0;j<numpnts(perfectSeq);j+=1)
					hold=selectnumber(allSeqAllAn[p][0]==i && allSeqAllAn[p][2]==perfectSeq[j],hold[p],allSeqAllAn[p][1])
				endfor
				sort hold,hold
				wavetransform zapNaNs,hold
				if(numpnts(hold)>=thr)
					trialToPerfect[i]=hold[thr-1]
				else
					trialToPerfect[i]=inf
				endif
			endif
		endfor
	else
		make /o/n=(dimsize(allAnContrib,0)) trialToPerfect=NaN
	endif
end

function didntReachPerfect(reOrderWv,which)
	wave reOrderWv
	variable which
	wave wv=$"allTrialToPerfect"+num2str(which)
	variable cntr,tot
	wavestats /q/rmd=[][1] reOrderWv
	variable howMany=v_max+1
	make /o/n=(howMany*2,dimsize(wv,1)) fracMissed=0
	variable i,j
	for(i=0;i<dimsize(fracMissed,1);i+=1)
		for(j=0;j<dimsize(reOrderWv,0);j+=1)
			if(wv[reOrderWv[j][0]][i]>=0)
				fracMissed[reOrderWv[j][1]][i]+=wv[reOrderWv[j][0]][i]==inf
				fracMissed[reOrderWv[j][1]+howMany][i]+=1
			endif
		endfor
	endfor
end

function doResampleAndFindK(reOrderWvList,full)
	string reOrderWvList
	variable full
	
	wave reOrderWv=$stringFromList(0,reOrderWvList)
	
	variable alph=0.5
	variable gam=0.05
	
	variable kVal=1
	
	make /o/n=(0) rejects,kVals
	
	if(getAllMatrixRows(reOrderWv))
	
		variable i,j
		do
			getRWpVals(kVal,alph)
			wave rwPvals
			duplicate /o rwPvals numP
			numP=rwPvals<alph
			make /o/n=(i+1) rejects,kVals
			rejects[i]=sum(numP)
			kVals[i]=kVal
			
			if(i==0)
				for(j=0;j<itemsInList(reOrderWvList);j+=1)
					duplicate /o/r=[numpnts(rwPvals)/itemsInList(reOrderWvList)*j,numpnts(rwPvals)/itemsInList(reOrderWvList)*(j+1)-1][] rwPvals $"rwPvalsFWER_"+stringFromList(j,reOrderWvList)
				endfor
			endif
			
			if(((kVal/gam)-1)==sum(numP))
				kVal+=1
			else
				kVal=ceil((sum(numP)+1)*gam)
			endif
			
			putRowsBackTogther()
			doUpdate
			i+=1
		while(sum(numP)>=(kVals[i-1]/gam-1))
		
		if(full)
			getRWpVals(kVals[i-1],1)
		endif
		
		for(j=0;j<itemsInList(reOrderWvList);j+=1)
			duplicate /o/r=[numpnts(rwPvals)/itemsInList(reOrderWvList)*j,numpnts(rwPvals)/itemsInList(reOrderWvList)*(j+1)-1][] rwPvals $"rwPvals_"+stringFromList(j,reOrderWvList)
		endfor
		
		killwaves /z numP
		
		killdatafolder /z resampDF
		getEveryCor(reOrderWv)
		saveExperiment
		return 1
	else
		print "too many comparisons, index matrix too small"
		return 0
	endif
end

function getRWpVals(k,maxPval)
	variable k,maxPval
	
	variable timer=startMStimer
	
	NVAR totalComp
	
	make /o/n=(totalComp) rwPvals=NaN
	variable lastP=-1
	variable i
	for(i=0;i<numpnts(rwPvals);i+=1)
		lastP=getOneRWpValue(rwPvals,k,lastP)
		if(lastP>=maxPval)
			rwPvals=selectnumber(rwPvals[p]^2>=0,lastP,rwPvals[p])
			break
		endif
		if(mod(i+1,1000)==0)
			printf "------Calculated %g p values------\r",i+1
		endif
	endfor
	
	print "------------------------------------"
	printf "        took %g s\r ",stopMStimer(timer)/1e6
	print "------------------------------------"
end

function putRowsBackTogther()
	
	variable timer=startMStimer
	
	print "--------------------------------------------"
	print "------puting row vectors back together------"
	
	wave wv=root:resampDF:data
	wave wvDone=root:resampDF:dataDone
	wave wvInd=root:resampDF:dataInd
	wave wvIndDone=root:resampDF:dataIndDone
	variable val1=numpnts(wv)
	variable val2=val1+numpnts(wvDone)
	
	make /o/n=(val2) root:resampDF:data
	make /o/i/u/n=(val2) root:resampDF:dataInd
	wv[val1,val2-1]=wvDone[p-val1]
	wvInd[val1,val2-1]=wvIndDone[p-val1]
	killwaves /z wvDone,wvIndDone
	sort /r wv,wv,wvInd
	
	variable i
	for(i=0;i<ResampleNumber;i+=1)
		wave wv=$"root:resampDF:row"+num2istr(i)
		wave wvDone=$"root:resampDF:rowDone"+num2istr(i)
		wave wvInd=$"root:resampDF:rowInd"+num2istr(i)
		wave wvIndDone=$"root:resampDF:rowIndDone"+num2istr(i)
		if(waveExists(wvDone))
			val1=numpnts(wv)
			val2=val1+numpnts(wvDone)
			make /o/n=(val2) $"root:resampDF:row"+num2istr(i)
			make /o/i/u/n=(val2) $"root:resampDF:rowInd"+num2istr(i)
			wv[val1,val2-1]=wvDone[p-val1]
			wvInd[val1,val2-1]=wvIndDone[p-val1]
			killwaves /z wvDone,wvIndDone
		endif
	endfor
	
	doParSort()
	
	printf "            took %g s\r ",stopMStimer(timer)/1e6
	print "---------------------------------------------"
	
	killwaves /z holdData,holdDataInd,holdDone,holdInd,holdIndDone,holdRow,holdRowInd
end

Function doParSort()
		
	Variable ncol=ResampleNumber
	Variable i,col,nthreads= ThreadProcessorCount
	Variable threadGroupID= ThreadGroupCreate(nthreads)
	
	for(col=0; col<ncol;)
		for(i=0; i<nthreads; i+=1)
			wave wv=$"root:resampDF:row"+num2istr(col)
			wave wvInd=$"root:resampDF:rowInd"+num2istr(col)
			ThreadStart threadGroupID,i,parSort(wv,wvInd)
			col+=1
			if( col>=ncol )
				break
			endif
		endfor
		
		do
			Variable threadGroupStatus = ThreadGroupWait(threadGroupID,100)
		while( threadGroupStatus != 0 )
	endfor
	Variable dummy= ThreadGroupRelease(threadGroupID)
End

ThreadSafe Function parSort(wv,wvInd)
	WAVE wv,wvInd
	
	sort /r wv,wv,wvInd
	
	return stopMSTimer(-2)		// time when we finished
End

function getOneRWpValue(pValWv,k,lastP)
	wave pValWv
	variable k,lastP
	make /o/n=(ResampleNumber) null=NaN,removeVals=NaN
	
	wave wv=$"root:resampDF:data"
	variable val=wv[0]
	wave wv=$"root:resampDF:dataInd"
	variable idx=wv[0]
	
	Variable ncol=ResampleNumber
	Variable i,col,nthreads= ThreadProcessorCount
	Variable threadGroupID= ThreadGroupCreate(nthreads)
	
	for(col=0; col<ncol;)
		for(i=0; i<nthreads; i+=1)
			wave rowWv=$"root:resampDF:row"+num2istr(col)
			wave idxWv=$"root:resampDF:rowInd"+num2istr(col)
			ThreadStart threadGroupID,i,nullAndIndex(rowWv,idxWv,null,removeVals,k,idx,col)
			col+=1
			if( col>=ncol )
				break
			endif
		endfor
		
		do
			Variable threadGroupStatus = ThreadGroupWait(threadGroupID,100)
		while( threadGroupStatus != 0 )
	endfor
	Variable dummy= ThreadGroupRelease(threadGroupID)

	null=(null>=val)
	variable pVal=(sum(null)+1)/(ResampleNumber+1)
	pValWv[idx]=max(pVal,lastP)
	
	removeEntries(removeVals)
	
	killwaves /z null,removeVals
	
	return pVal
end

threadsafe function nullAndIndex(rowWv,idxWv,nullWv,rmvWv,kVal,idx,col)
	wave rowWv,idxWv,nullWv,rmvWv
	variable kVal,idx,col
	
	nullWv[col]=rowWv[kVal-1]
	
	findvalue /i=(idx) idxWv
	rmvWv[col]=v_value
	
	return stopMSTimer(-2)		// time when we finished
end

function removeEntries(rmvWv)
	wave rmvWv
	
	wave wv=$"root:resampDF:data"
	make /o/n=(1) holdVal=wv[0]
	deletepoints 0,1,wv
	wave wv=$"root:resampDF:dataInd"
	make /o/u/i/n=(1) holdInd=wv[0]
	deletepoints 0,1,wv
	
	concatenate /NP=0 "holdVal;",root:resampDF:dataDone
	concatenate /NP=0 "holdInd;",root:resampDF:dataIndDone
	
	variable i
	for(i=0;i<numpnts(rmvWv);i+=1)
		if(rmvWv[i]>=0)
			wave wv=$"root:resampDF:row"+num2istr(i)
			make /o/n=(1) holdVal=wv[rmvWv[i]]
			deletepoints rmvWv[i],1,wv
			wave wv=$"root:resampDF:rowInd"+num2istr(i)
			make /o/u/i/n=(1) holdInd=wv[rmvWv[i]]
			deletepoints rmvWv[i],1,wv
			
			concatenate /NP=0 "holdVal;",$"root:resampDF:rowDone"+num2istr(i)
			concatenate /NP=0 "holdInd;",$"root:resampDF:rowIndDone"+num2istr(i)
		endif
	endfor
	
	killwaves /z holdVal,holdInd
end

function getAllMatrixRows(reOrderWv)
	wave reOrderWv
	
	variable i,j,k
	
	variable conting,len	
	make /o/n=(itemsInList(ContingList)*SeqLenMax,3) seqCntInfo=NaN
	for(i=0;i<itemsInList(ContingList);i+=1)
		conting=str2num(stringFromList(i,ContingList))
		for(j=0;j<SeqLenMax;j+=1)
			len=j+1
			wave cntWv=$"root:seqDF:seqCnts_"+num2str(conting)+"_"+num2str(len)
			seqCntInfo[k][0]=conting
			seqCntInfo[k][1]=len
			seqCntInfo[k][2]=dimsize(cntWv,1)
			k+=1
		endfor
	endfor
	
	wavestats /q/rmd=[][2] seqCntInfo
	variable totSeq=v_sum
	variable /g totalComp=totSeq*2
	
	make /o/n=(dimsize(reOrderWv,0)) subNums=reOrderWv[p][0],covNums=reOrderWv[p][1]
	
	if((totalComp)>=2^32)
		print "ERROR!!!! Too Many comparisons"
		return 0
	endif
	
	killdatafolder /z resampDF
	newdataFolder resampDF
	
	getMatrixRow(seqCntInfo,subNums,covNums)
	wave oneRow,rowInd
	duplicate /o oneRow $"root:resampDF:data"
	duplicate /o rowInd $"root:resampDF:dataInd"
	
	printf "--------Data Studentization done---: %g sequences\r",totSeq
	
	make /o/n=(numpnts(subNums)) rand,order
	for(i=0;i<ResampleNumber;i+=1)
		setRandomSeed /BETR (i+1)/ResampleNumber
		order=subNums[p]
		rand=gnoise(1)
		sort rand,order
		getMatrixRow(seqCntInfo,order,covNums)
		duplicate /o oneRow $"root:resampDF:row"+num2istr(i)
		duplicate /o rowInd $"root:resampDF:rowInd"+num2istr(i)
		
		if(mod(i+1,1000)==0)
			printf "--------%g resamples done----------: %g sequences\r",i+1,totSeq
		endif
	endfor
	
	killwaves /z oneRow,rowInd,rand,order,subNums,covNums
			
	return 1
end

//assumes fewer than 2^32 comparisons
Function getMatrixRow(seqSizeWv,subWv,covWv)
	wave seqSizeWv,subWv,covWv
	
	wave allAnContrib
	
	Variable j
	Variable numSeq,seqNum
	Variable i,nthreads= ThreadProcessorCount
	Variable threadGroupID= ThreadGroupCreate(nthreads)
	variable cnt=0
	
	NVAR totalComp
	make /o/n=(totalComp) oneRow=NaN
	
	for(j=0;j<dimsize(seqSizeWv,0);j+=1)
		numSeq=seqSizeWv[j][2]
		wave cntWv=$"root:seqDF:seqCnts_"+num2str(seqSizeWv[j][0])+"_"+num2str(seqSizeWv[j][1])
		make /o/n=(numpnts(subWv)) pnts=selectnumber(allAnContrib[subWv[p]][seqSizeWv[j][0]]==1,NaN,p)
		wavetransform zapNaNs,pnts
		make /o/n=(numpnts(pnts)) subWv4Conting=subWv[pnts[p]],covWv4Conting=covWv[pnts[p]]
		for(seqNum=0; seqNum<numSeq;)
			for(i=0; i<nthreads; i+=1)
				ThreadStart threadGroupID,i,getSingleMatrixEntry(cntWv,subWv4Conting,covWv4Conting,oneRow,seqNum,cnt+seqNum)
				seqNum+=1
				if( seqNum>=numSeq )
					break
				endif
			endfor
			
			do
				Variable threadGroupStatus = ThreadGroupWait(threadGroupID,100)
			while( threadGroupStatus != 0 )
		endfor
		cnt+=numSeq
	endfor
	Variable dummy= ThreadGroupRelease(threadGroupID)
	
	make /o/u/i/n=(numpnts(oneRow)) rowInd
	
	rowInd=p
	sort /r oneRow,oneRow,rowInd
	wavestats /q oneRow
	deletepoints 0,v_numNaNs,oneRow,rowInd
End

ThreadSafe Function getSingleMatrixEntry(wv,subWv,covWv,outWv,which,cnt)
	WAVE wv,subWv,covWv,outWv
	Variable which,cnt
	
	make /o/n=(numpnts(subWv)) cntVals=wv[subWv[p]][which]
	
	variable top=0,bot1=0,bot2=0
	wavestats /q/m=1 cntVals
	variable mn1=v_avg
	wavestats /q/m=1 covWv
	variable mn2=v_avg
	variable i
	for(i=0;i<numpnts(cntVals);i+=1)
		top+=(cntVals[i]-mn1)^2*(covWv[i]-mn2)^2
		bot1+=(cntVals[i]-mn1)^2
		bot2+=(covWv[i]-mn2)^2
	endfor
	variable val=statscorrelation(cntVals,covWv)
	if(val^2>=0)
		if(val>0)
			outWv[cnt*2+0]=val/(sqrt(top)/(sqrt(bot1)*sqrt(bot2)))
			outWv[cnt*2+1]=NaN
		elseif(val<0)
			outWv[cnt*2+0]=NaN
			outWv[cnt*2+1]=-val/(sqrt(top)/(sqrt(bot1)*sqrt(bot2)))
		else
			outWv[cnt*2+0]=NaN
			outWv[cnt*2+1]=NaN
		endif
	else
		outWv[cnt*2+0]=NaN
		outWv[cnt*2+1]=NaN
	endif
		
	return stopMSTimer(-2)		// time when we finished
End

function getEveryCor(reOrderWv)
	wave reOrderWv
	wave seqCntInfo
	
	wavestats /q/rmd=[][1] reOrderWv
	variable mx=v_max+1
	
	wavestats /q/rmd=[][2] seqCntInfo
	make /o/n=(v_sum,5) everyCor=NaN
	
	Variable j
	Variable numSeq,seqNum
	Variable i,nthreads= ThreadProcessorCount
	Variable threadGroupID= ThreadGroupCreate(nthreads)
	variable cnt=0
		
	for(j=0;j<dimsize(seqCntInfo,0);j+=1)
		numSeq=seqCntInfo[j][2]
		wave cntWv=$"root:seqDF:seqCnts_"+num2str(seqCntInfo[j][0])+"_"+num2str(seqCntInfo[j][1])
		everyCor[cnt,cnt+numSeq-1][0,1]=seqCntInfo[j][q]
		everyCor[cnt,cnt+numSeq-1][2]=p-cnt
		for(seqNum=0; seqNum<numSeq;)
			for(i=0; i<nthreads; i+=1)
				ThreadStart threadGroupID,i,getCors(cntWv,reOrderWv,everyCor,seqNum,cnt+seqNum)
				seqNum+=1
				if( seqNum>=numSeq )
					break
				endif
			endfor
			
			do
				Variable threadGroupStatus = ThreadGroupWait(threadGroupID,100)
			while( threadGroupStatus != 0 )
		endfor
		cnt+=numSeq
	endfor
	Variable dummy= ThreadGroupRelease(threadGroupID)
	
	duplicate /o everyCor $"everyCor_"+nameOfWave(reOrderWv)
end

ThreadSafe Function getCors(cntWv,rWv,cWv,which,cnt)
	WAVE cntWv,rWv,cWv
	Variable which,cnt
	
	make /o/n=(dimsize(rWv,0)) pnts=rWv[p][0],cov=rWv[p][1]
	wavetransform zapNaNs,pnts
	if(numpnts(pnts)>0)
		make /o/n=(numpnts(pnts)) hold=cntWv[pnts[p]][which]
		cWv[cnt][3]=statscorrelation(hold,cov)
		wavestats /q/m=1 hold
		cWv[cnt][4]=v_avg
	else
		cWv[cnt][3]=NaN
		cWv[cnt][4]=NaN
	endif
		
	return stopMSTimer(-2)		// time when we finished
End