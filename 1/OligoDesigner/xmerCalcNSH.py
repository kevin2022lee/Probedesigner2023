from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
class CalcNSH:
    def __init__(self):
        pass
    '''This a Class that can calculate the NSH and LCS'''
    def lcs(self,first,second):  
        first_length = len(first)  
        second_length = len(second)  
        size = 0  
        x = 0  
        y = 0  
        matrix = [range(second_length) for x in range(first_length)]  
        #print matrix  
        for i in range(first_length):  
            for j in range(second_length):  
                #print i,j  
                if first[i] == second[j]:  
                    if i - 1 >= 0 and j - 1 >=0:  
                        matrix[i][j] = matrix[i-1][j-1] + 1  
                    else:  
                        matrix[i][j] = 1  
                    if matrix[i][j] > size:  
                        size = matrix[i][j]  
                        x = j  
                        y = i  
                else:  
                    matrix[i][j] = 0  
        #print matrix  
        #print size,x,y   
  
        return second[x-size+1:x+1] 
    ############################################
    def QuantimatxmerCalc(self,Calc_Seq):
#########################QuantiMAT2.0 universal file#####################################################    
        uni_Aleader_seq=Seq("AAAACGGTAACTTCTTTATGCTTTGACTCAG", IUPAC.unambiguous_dna)
        uni_Aarms_seq=Seq("ATCTCAGTCTCGTTAATGGATTCCT", IUPAC.unambiguous_dna)
        uni_AP_seq=Seq("AAGTACGACAACCACATCTT", IUPAC.unambiguous_dna)
        uni_PSCP_seq=Seq("CACTTCACTTTCTTTCCAAGAG", IUPAC.unambiguous_dna)
        
        ########LCS#############
        Calc_Seq=Seq(Calc_Seq, IUPAC.unambiguous_dna).upper()
        x4merlcs_Aleader=self.lcs(str(uni_Aleader_seq),str(Calc_Seq))
        x4merlcs_Aarms=self.lcs(str(uni_Aarms_seq),str(Calc_Seq))
        x4merlcs_AP=self.lcs(str(uni_AP_seq),str(Calc_Seq))
        x4merlcs_PSCP=self.lcs(str(uni_PSCP_seq),str(Calc_Seq))
###########Server as CE Probe Weighting Factor##########################
        WF_CEtoLeaders=3
        WF_CEtoAMParms=20
        WF_CEtoAP=1
        WF_CEtoPSCP=0
###########Server as LE Probe Weighting Factor##########################
        WF_LEtoLeaders=0
        WF_LEtoAMParms=0
        WF_LEtoAP=0
        WF_LEtoPSCP=12
############################################################################
        if len(x4merlcs_Aleader) >=4:
            score_x4mer_Aleader=[]
            data_Aleader={}
            i=0
            while(i<len(x4merlcs_Aleader)-3):
                score_x4mer_Aleader.append(self.x4merScore(x4merlcs_Aleader[i:i+4]))
                i=i+1
            data_Aleader['score_x4mer_Aleader']=score_x4mer_Aleader
            NSH_Score_Aleader_SACE=sum(data_Aleader['score_x4mer_Aleader'])*WF_CEtoLeaders
            NSH_Score_Aleader_SALE=sum(data_Aleader['score_x4mer_Aleader'])*WF_LEtoLeaders
        else:
            NSH_Score_Aleader_SACE=0
            NSH_Score_Aleader_SALE=0
###############Aarms xmer Score############################################
        if len(x4merlcs_Aarms) >=4:
            i=0
            score_x4mer_Aarms=[]
            data_Aarms={}
            while(i<len(x4merlcs_Aarms)-3):
                score_x4mer_Aarms.append(self.x4merScore(x4merlcs_Aarms[i:i+4]))
                i=i+1
            data_Aarms['score_x4mer_Aarms']=score_x4mer_Aarms
            NSH_Score_Aarms_SACE=sum(data_Aarms['score_x4mer_Aarms'])*WF_CEtoAMParms
            NSH_Score_Aarms_SALE=sum(data_Aarms['score_x4mer_Aarms'])*WF_LEtoAMParms
        else:
            NSH_Score_Aarms_SACE=0
            NSH_Score_Aarms_SALE=0
###############APs xmer Score##############################################
        if len(x4merlcs_AP) >=4:
            i=0
            score_x4mer_AP=[]
            data_AP={}
            while(i<len(x4merlcs_AP)-3):
                score_x4mer_AP.append(self.x4merScore(x4merlcs_AP[i:i+4]))
                i=i+1
            data_AP['score_x4mer_AP']=score_x4mer_AP
            NSH_Score_AP_SACE=sum(data_AP['score_x4mer_AP'])*WF_CEtoAP
            NSH_Score_AP_SALE=sum(data_AP['score_x4mer_AP'])*WF_LEtoAP
        else:
            NSH_Score_AP_SACE=0
            NSH_Score_AP_SALE=0
##############PSCP x-mer Score###############################
        if len(x4merlcs_PSCP) >=4:
            i=0
            score_x4mer_PSCP=[]
            data_PSCP={}
            while(i<len(x4merlcs_PSCP)-3):
                score_x4mer_PSCP.append(self.x4merScore(x4merlcs_PSCP[i:i+4]))
                i=i+1
            data_PSCP['score_x4mer_PSCP']=score_x4mer_PSCP
            NSH_Score_PSCP_SACE=sum(data_PSCP['score_x4mer_PSCP'])*WF_CEtoPSCP
            NSH_Score_PSCP_SALE=sum(data_PSCP['score_x4mer_PSCP'])*WF_LEtoPSCP
        else:
            NSH_Score_PSCP_SACE=0
            NSH_Score_PSCP_SALE=0
#######END############
        Total_NSH_SACE= NSH_Score_Aleader_SACE+NSH_Score_Aarms_SACE+NSH_Score_AP_SACE+NSH_Score_PSCP_SACE
        Total_NSH_SALE= NSH_Score_Aleader_SALE+NSH_Score_Aarms_SALE+NSH_Score_AP_SALE+NSH_Score_PSCP_SALE
        return [Total_NSH_SACE,Total_NSH_SALE]
    
    def QuantiplexxmerCalc(self,Calc_Seq):
#########################Quantiplex 2.0 universal file#####################################################    
        uni_Aleader_seq=Seq("CGGCATAGCAGCGCGCATACTCTTCCGGTCGCCCAATGGTCC", IUPAC.unambiguous_dna)
        uni_Aarms_seq=Seq("CGTTGTCCCTAGGGCCGTGGATTTTT", IUPAC.unambiguous_dna)
        uni_PSCP_seq=Seq("CACTTCACTTTCTTTCCAAGAG", IUPAC.unambiguous_dna)
        uni_AP_seq=Seq("GCACTTGGTACGGCGCTGACTTTTTTT", IUPAC.unambiguous_dna)
        
        ########LCS#############
        Calc_Seq=Seq(Calc_Seq, IUPAC.unambiguous_dna).upper()
        x4merlcs_Aleader=self.lcs(str(uni_Aleader_seq),str(Calc_Seq))
        x4merlcs_Aarms=self.lcs(str(uni_Aarms_seq),str(Calc_Seq))
        x4merlcs_AP=self.lcs(str(uni_AP_seq),str(Calc_Seq))
        x4merlcs_PSCP=self.lcs(str(uni_PSCP_seq),str(Calc_Seq))
###########Server as CE Probe Weighting Factor##########################
        WF_CEtoLeaders=3
        WF_CEtoAMParms=14
        WF_CEtoAP=1
        WF_CEtoPSCP=0
###########Server as LE Probe Weighting Factor##########################
        WF_LEtoLeaders=0
        WF_LEtoAMParms=0
        WF_LEtoAP=0
        WF_LEtoPSCP=6
############################################################################
        if len(x4merlcs_Aleader) >=4:
            score_x4mer_Aleader=[]
            data_Aleader={}
            i=0
            while(i<len(x4merlcs_Aleader)-3):
                score_x4mer_Aleader.append(self.x4merScore(x4merlcs_Aleader[i:i+4]))
                i=i+1
            data_Aleader['score_x4mer_Aleader']=score_x4mer_Aleader
            NSH_Score_Aleader_SACE=sum(data_Aleader['score_x4mer_Aleader'])*WF_CEtoLeaders
            NSH_Score_Aleader_SALE=sum(data_Aleader['score_x4mer_Aleader'])*WF_LEtoLeaders
        else:
            NSH_Score_Aleader_SACE=0
            NSH_Score_Aleader_SALE=0
###############Aarms xmer Score############################################
        if len(x4merlcs_Aarms) >=4:
            i=0
            score_x4mer_Aarms=[]
            data_Aarms={}
            while(i<len(x4merlcs_Aarms)-3):
                score_x4mer_Aarms.append(self.x4merScore(x4merlcs_Aarms[i:i+4]))
                i=i+1
            data_Aarms['score_x4mer_Aarms']=score_x4mer_Aarms
            NSH_Score_Aarms_SACE=sum(data_Aarms['score_x4mer_Aarms'])*WF_CEtoAMParms
            NSH_Score_Aarms_SALE=sum(data_Aarms['score_x4mer_Aarms'])*WF_LEtoAMParms
        else:
            NSH_Score_Aarms_SACE=0
            NSH_Score_Aarms_SALE=0
###############APs xmer Score##############################################
        if len(x4merlcs_AP) >=4:
            i=0
            score_x4mer_AP=[]
            data_AP={}
            while(i<len(x4merlcs_AP)-3):
                score_x4mer_AP.append(self.x4merScore(x4merlcs_AP[i:i+4]))
                i=i+1
            data_AP['score_x4mer_AP']=score_x4mer_AP
            NSH_Score_AP_SACE=sum(data_AP['score_x4mer_AP'])*WF_CEtoAP
            NSH_Score_AP_SALE=sum(data_AP['score_x4mer_AP'])*WF_LEtoAP
        else:
            NSH_Score_AP_SACE=0
            NSH_Score_AP_SALE=0
##############PSCP x-mer Score###############################
        if len(x4merlcs_PSCP) >=4:
            i=0
            score_x4mer_PSCP=[]
            data_PSCP={}
            while(i<len(x4merlcs_PSCP)-3):
                score_x4mer_PSCP.append(self.x4merScore(x4merlcs_PSCP[i:i+4]))
                i=i+1
            data_PSCP['score_x4mer_PSCP']=score_x4mer_PSCP
            NSH_Score_PSCP_SACE=sum(data_PSCP['score_x4mer_PSCP'])*WF_CEtoPSCP
            NSH_Score_PSCP_SALE=sum(data_PSCP['score_x4mer_PSCP'])*WF_LEtoPSCP
        else:
            NSH_Score_PSCP_SACE=0
            NSH_Score_PSCP_SALE=0
#######END############
        Total_NSH_SACE= NSH_Score_Aleader_SACE+NSH_Score_Aarms_SACE+NSH_Score_AP_SACE+NSH_Score_PSCP_SACE
        Total_NSH_SALE= NSH_Score_Aleader_SALE+NSH_Score_Aarms_SALE+NSH_Score_AP_SALE+NSH_Score_PSCP_SALE
        return [Total_NSH_SACE,Total_NSH_SALE]
    
    def Quantimat3xmerCalc(self,Calc_Seq):
#########################QuantiMAT 3.0 universal file#####################################################    
        
        uni_AB_seq=Seq("ACTGTCAAATCT", IUPAC.unambiguous_dna)
        uni_CD_seq=Seq("TTCAAATGCTCA", IUPAC.unambiguous_dna)
        uni_EF_seq=Seq("TCACGAATCTAC", IUPAC.unambiguous_dna)
        uni_GH_seq=Seq("GTGCTTGTGCCT", IUPAC.unambiguous_dna)
        uni_IJ_seq=Seq("GTCATAGCAAAT", IUPAC.unambiguous_dna)
        uni_KL_seq=Seq("ACTTATCGCATC", IUPAC.unambiguous_dna)
        uni_MN_seq=Seq("CCTGTTTATTGC", IUPAC.unambiguous_dna)
        uni_PSCP_seq=Seq("CACTTCACTTTCTTTCCAAGAG", IUPAC.unambiguous_dna)
        
        uni_CE_TAIL_seq=Seq("CTCTTGGAAAGAAAGTGAAGTG", IUPAC.unambiguous_dna)
        uni_LE_TAIL_seq=Seq("JTATJCGCJCTGFTATJCCGT", IUPAC.unambiguous_dna)
        
        ########LCS############# 3.0 NSH all from Linkers   2023.5.10##########fuck code#############
        Calc_Seq=Seq(Calc_Seq, IUPAC.unambiguous_dna).upper()
        
        x4merlcs_AB=self.lcs(str(uni_AB_seq),str(Calc_Seq))
        x4merlcs_CD=self.lcs(str(uni_CD_seq),str(Calc_Seq))
        x4merlcs_EF=self.lcs(str(uni_EF_seq),str(Calc_Seq))
        x4merlcs_GH=self.lcs(str(uni_GH_seq),str(Calc_Seq))
        x4merlcs_IJ=self.lcs(str(uni_IJ_seq),str(Calc_Seq))
        x4merlcs_KL=self.lcs(str(uni_KL_seq),str(Calc_Seq))
        x4merlcs_MN=self.lcs(str(uni_MN_seq),str(Calc_Seq))
        x4merlcs_PSCP=self.lcs(str(uni_PSCP_seq),str(Calc_Seq))
###########Server as CE Probe Weighting Factor##########################
        WF_CEtoAB=8
        WF_CEtoCD=8
        WF_CEtoEF=8
        WF_CEtoGH=8
        WF_CEtoIJ=8
        WF_CEtoKL=8
        WF_CEtoMN=8
        WF_CEtoPSCP=0
###########Server as LE Probe Weighting Factor##########################
        WF_LEtoAB=0
        WF_LEtoCD=0
        WF_LEtoEF=0
        WF_LEtoGH=0
        WF_LEtoIJ=0
        WF_LEtoKL=0
        WF_LEtoMN=0
        WF_LEtoPSCP=8
############################################################################
        if len(x4merlcs_AB) >=4:
            score_x4mer_AB=[]
            data_AB={}
            i=0
            while(i<len(x4merlcs_AB)-3):
                score_x4mer_AB.append(self.x4merScore(x4merlcs_AB[i:i+4]))
                i=i+1
            data_AB['score_x4mer_AB']=score_x4mer_AB
            NSH_Score_AB_SACE=sum(data_AB['score_x4mer_AB'])*WF_CEtoAB
            NSH_Score_AB_SALE=sum(data_AB['score_x4mer_AB'])*WF_LEtoAB
        else:
            NSH_Score_AB_SACE=0
            NSH_Score_AB_SALE=0
###############LinkerAB xmer Score############################################
        if len(x4merlcs_CD) >=4:
            score_x4mer_CD=[]
            data_CD={}
            i=0
            while(i<len(x4merlcs_CD)-3):
                score_x4mer_CD.append(self.x4merScore(x4merlcs_CD[i:i+4]))
                i=i+1
            data_CD['score_x4mer_CD']=score_x4mer_CD
            NSH_Score_CD_SACE=sum(data_CD['score_x4mer_CD'])*WF_CEtoCD
            NSH_Score_CD_SALE=sum(data_AB['score_x4mer_CD'])*WF_LEtoCD
        else:
            NSH_Score_CD_SACE=0
            NSH_Score_CD_SALE=0
###############LinkerCD xmer Score##############################################
        if len(x4merlcs_EF) >=4:
            score_x4mer_EF=[]
            data_EF={}
            i=0
            while(i<len(x4merlcs_EF)-3):
                score_x4mer_EF.append(self.x4merScore(x4merlcs_EF[i:i+4]))
                i=i+1
            data_EF['score_x4mer_EF']=score_x4mer_EF
            NSH_Score_EF_SACE=sum(data_EF['score_x4mer_EF'])*WF_CEtoEF
            NSH_Score_EF_SALE=sum(data_EF['score_x4mer_EF'])*WF_LEtoEF
        else:
            NSH_Score_EF_SACE=0
            NSH_Score_EF_SALE=0
##############LinkerEF x-mer Score###############################
        if len(x4merlcs_GH) >=4:
            score_x4mer_GH=[]
            data_GH={}
            i=0
            while(i<len(x4merlcs_GH)-3):
                score_x4mer_GH.append(self.x4merScore(x4merlcs_GH[i:i+4]))
                i=i+1
            data_GH['score_x4mer_GH']=score_x4mer_GH
            NSH_Score_GH_SACE=sum(data_GH['score_x4mer_GH'])*WF_CEtoGH
            NSH_Score_GH_SALE=sum(data_GH['score_x4mer_GH'])*WF_LEtoGH
        else:
            NSH_Score_GH_SACE=0
            NSH_Score_GH_SALE=0
##############LinkerGH x-mer Score###############################
        if len(x4merlcs_IJ) >=4:
            score_x4mer_IJ=[]
            data_IJ={}
            i=0
            while(i<len(x4merlcs_IJ)-3):
                score_x4mer_IJ.append(self.x4merScore(x4merlcs_IJ[i:i+4]))
                i=i+1
            data_IJ['score_x4mer_IJ']=score_x4mer_IJ
            NSH_Score_IJ_SACE=sum(data_IJ['score_x4mer_IJ'])*WF_CEtoIJ
            NSH_Score_IJ_SALE=sum(data_IJ['score_x4mer_IJ'])*WF_LEtoIJ
        else:
            NSH_Score_IJ_SACE=0
            NSH_Score_IJ_SALE=0
##############LinkerIJ x-mer Score###############################
        if len(x4merlcs_KL) >=4:
            score_x4mer_KL=[]
            data_KL={}
            i=0
            while(i<len(x4merlcs_KL)-3):
                score_x4mer_KL.append(self.x4merScore(x4merlcs_KL[i:i+4]))
                i=i+1
            data_KL['score_x4mer_KL']=score_x4mer_KL
            NSH_Score_KL_SACE=sum(data_KL['score_x4mer_KL'])*WF_CEtoKL
            NSH_Score_KL_SALE=sum(data_KL['score_x4mer_KL'])*WF_LEtoKL
        else:
            NSH_Score_KL_SACE=0
            NSH_Score_KL_SALE=0
##############LinkerKL x-mer Score###############################
        if len(x4merlcs_MN) >=4:
            score_x4mer_MN=[]
            data_MN={}
            i=0
            while(i<len(x4merlcs_MN)-3):
                score_x4mer_MN.append(self.x4merScore(x4merlcs_MN[i:i+4]))
                i=i+1
            data_MN['score_x4mer_MN']=score_x4mer_MN
            NSH_Score_MN_SACE=sum(data_MN['score_x4mer_MN'])*WF_CEtoMN
            NSH_Score_MN_SALE=sum(data_MN['score_x4mer_MN'])*WF_LEtoMN
        else:
            NSH_Score_MN_SACE=0
            NSH_Score_MN_SALE=0
##############LinkerMN x-mer Score###############################
        if len(x4merlcs_PSCP) >=4:
            i=0
            score_x4mer_PSCP=[]
            data_PSCP={}
            while(i<len(x4merlcs_PSCP)-3):
                score_x4mer_PSCP.append(self.x4merScore(x4merlcs_PSCP[i:i+4]))
                i=i+1
            data_PSCP['score_x4mer_PSCP']=score_x4mer_PSCP
            NSH_Score_PSCP_SACE=sum(data_PSCP['score_x4mer_PSCP'])*WF_CEtoPSCP
            NSH_Score_PSCP_SALE=sum(data_PSCP['score_x4mer_PSCP'])*WF_LEtoPSCP
        else:
            NSH_Score_PSCP_SACE=0
            NSH_Score_PSCP_SALE=0
#######END###################################################################
        Total_NSH_SACE= NSH_Score_AB_SACE+NSH_Score_CD_SACE+NSH_Score_EF_SACE+NSH_Score_GH_SACE+NSH_Score_IJ_SACE+NSH_Score_KL_SACE+NSH_Score_MN_SACE+NSH_Score_PSCP_SACE
        Total_NSH_SALE= NSH_Score_AB_SALE+NSH_Score_CD_SALE+NSH_Score_EF_SALE+NSH_Score_GH_SALE+NSH_Score_IJ_SALE+NSH_Score_KL_SALE+NSH_Score_MN_SALE+NSH_Score_PSCP_SALE
        return [Total_NSH_SACE,Total_NSH_SALE]

    def x4merScore(self,seq):
        SumAT=seq.count("A")+seq.count("T")
        SumGC=seq.count("G")+seq.count("C")
        Score=round((0.5*SumAT+1.0*SumGC)/4,3)
        return Score