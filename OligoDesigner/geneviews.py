﻿#coding:utf-8
from django.template import RequestContext
from django.shortcuts import render_to_response
from django.http import HttpResponse,HttpResponseRedirect,response
from django.core.paginator import Paginator
from django.views.decorators.csrf import csrf_protect,csrf_exempt
import time
from datetime import datetime
from models import *
from django.core.context_processors import request
import pymysql


local='www.probedesigner.cn'
thisyear=time.strftime('%Y',time.localtime(time.time()))


def genesearch(request):
    genes={} 
    if request.GET.get('genetype')=='Human':
        genes=GeneInfo.objects.all().reverse()[:100]  
    if request.GET.get('genetype')=='Mouse':
        genes=GeneInfo1.objects.all().reverse()[:100]
    if request.GET.get('genetype')=='Rat':
        genes=GeneInfo2.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='AtArabidopsis':
        genes=GeneInfo3.objects.all().reverse()[:100]   
    if request.GET.get('genetype')=='Celegans':
        genes=GeneInfo4.objects.all().reverse()[:100]  
    if request.GET.get('genetype')=='Fruitfly':
        genes=GeneInfo5.objects.all().reverse()[:100]
    if request.GET.get('genetype')=='Bovine':
        genes=GeneInfo6.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Dog':
        genes=GeneInfo7.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Chinesehamster':
        genes=GeneInfo8.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Goat':
        genes=GeneInfo9.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Guineapig':
        genes=GeneInfo10.objects.all().reverse()[:100] 
    #####################################################
    if request.GET.get('genetype')=='Zebrafish':
        genes=GeneInfo11.objects.all().reverse()[:100]  
    if request.GET.get('genetype')=='Horse':
        genes=GeneInfo12.objects.all().reverse()[:100]
    if request.GET.get('genetype')=='Chicken':
        genes=GeneInfo13.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Soybean':
        genes=GeneInfo14.objects.all().reverse()[:100]   
    if request.GET.get('genetype')=='Nakedmolerat':
        genes=GeneInfo15.objects.all().reverse()[:100]  
    if request.GET.get('genetype')=='CynomolgusMonkey':
        genes=GeneInfo16.objects.all().reverse()[:100]
    if request.GET.get('genetype')=='Sheep':
        genes=GeneInfo17.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Rabbit':
        genes=GeneInfo18.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Rice':
        genes=GeneInfo19.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Rhesusmonkeyhamster':
        genes=GeneInfo20.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Bakersyeast':
        genes=GeneInfo21.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Fissionyeast':
        genes=GeneInfo22.objects.all().reverse()[:100]  
    if request.GET.get('genetype')=='Pig':
        genes=GeneInfo23.objects.all().reverse()[:100]
    if request.GET.get('genetype')=='Breadwheat':
        genes=GeneInfo24.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Winegrape':
        genes=GeneInfo25.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Westernclawedfrog':
        genes=GeneInfo26.objects.all().reverse()[:100] 
    if request.GET.get('genetype')=='Maize':
        genes=GeneInfo27.objects.all().reverse()[:100] 


    stype=request.GET.get("genetype",'Human')
    return render_to_response('genedatabase/genesearch.html',{
                                    'local':local,
                                    'thisyear':thisyear,
                                    'genes':genes,
                                    'stype':stype
        },context_instance=RequestContext(request))

@csrf_protect
def gdbsearch(request):
    if request.method=="POST":
        genename=request.POST["search_text"] 
        if genename=="*" or genename=="":
            genes=GeneInfo.objects.filter(genename="GAPDH")
        else:
            table_num=request.POST["species"]
            if table_num=="":
                genes=GeneInfo.objects.filter(genename=genename)
            if  table_num=="1":
                genes=GeneInfo1.objects.filter(genename=genename)
            if  table_num=="2":
                genes=GeneInfo2.objects.filter(genename=genename)
            if  table_num=="3":
                genes=GeneInfo3.objects.filter(genename=genename)
            if  table_num=="4":
                genes=GeneInfo4.objects.filter(genename=genename)    
            if  table_num=="5":
                genes=GeneInfo5.objects.filter(genename=genename)
            if  table_num=="6":
                genes=GeneInfo6.objects.filter(genename=genename)
            if  table_num=="7":
                genes=GeneInfo7.objects.filter(genename=genename)
            if  table_num=="8":
                genes=GeneInfo8.objects.filter(genename=genename)
            if  table_num=="9":
                genes=GeneInfo9.objects.filter(genename=genename)
                ##########################
            if  table_num=="10":
                genes=GeneInfo10.objects.filter(genename=genename)
            if  table_num=="11":
                genes=GeneInfo11.objects.filter(genename=genename)
            if  table_num=="12":
                genes=GeneInfo12.objects.filter(genename=genename)
            if  table_num=="13":
                genes=GeneInfo13.objects.filter(genename=genename)    
            if  table_num=="14":
                genes=GeneInfo14.objects.filter(genename=genename)
            if  table_num=="15":
                genes=GeneInfo15.objects.filter(genename=genename)
            if  table_num=="16":
                genes=GeneInfo16.objects.filter(genename=genename)
            if  table_num=="17":
                genes=GeneInfo17.objects.filter(genename=genename)
            if  table_num=="18":
                genes=GeneInfo18.objects.filter(genename=genename)
                #######################################
            if  table_num=="19":
                genes=GeneInfo19.objects.filter(genename=genename)
            if  table_num=="20":
                genes=GeneInfo20.objects.filter(genename=genename)
            if  table_num=="21":
                genes=GeneInfo21.objects.filter(genename=genename)
            if  table_num=="22":
                genes=GeneInfo22.objects.filter(genename=genename)    
            if  table_num=="23":
                genes=GeneInfo23.objects.filter(genename=genename)
            if  table_num=="24":
                genes=GeneInfo24.objects.filter(genename=genename)
            if  table_num=="25":
                genes=GeneInfo25.objects.filter(genename=genename)
            if  table_num=="26":
                genes=GeneInfo26.objects.filter(genename=genename)
            if  table_num=="27":
                genes=GeneInfo27.objects.filter(genename=genename)
                
        return render_to_response('genedatabase/geneshow.html',{
                 'genes':genes,
                 'local':local,
                 'thisyear':thisyear
                 },context_instance=RequestContext(request))


def geneshow(request,specy,id):
    if specy=='Human':
        genes=GeneInfo.objects.filter(id__iexact=id)
    if specy=='Mouse':
        genes=GeneInfo1.objects.filter(id__iexact=id)
    if specy=='Rat':
        genes=GeneInfo2.objects.filter(id__iexact=id) 
    if specy=='AtArabidopsis':
        genes=GeneInfo3.objects.filter(id__iexact=id)
    if specy=='Celegans':
        genes=GeneInfo4.objects.filter(id__iexact=id)
    if specy=='Fruitfly':
        genes=GeneInfo5.objects.filter(id__iexact=id)
    if specy=='Bovine':
        genes=GeneInfo6.objects.filter(id__iexact=id) 
    if specy=='Dog':
        genes=GeneInfo7.objects.filter(id__iexact=id)
    if specy=='Chinesehamster':
        genes=GeneInfo8.objects.filter(id__iexact=id)
    if specy=='Goat':
        genes=GeneInfo9.objects.filter(id__iexact=id)
    ####################################################
    if specy=='Guineapig':
        genes=GeneInfo10.objects.filter(id__iexact=id)
    if specy=='Zebrafish':
        genes=GeneInfo11.objects.filter(id__iexact=id)
    if specy=='Horse':
        genes=GeneInfo12.objects.filter(id__iexact=id) 
    if specy=='Chicken':
        genes=GeneInfo13.objects.filter(id__iexact=id)
    if specy=='Soybean':
        genes=GeneInfo14.objects.filter(id__iexact=id)
    if specy=='Nakedmolerat':
        genes=GeneInfo15.objects.filter(id__iexact=id)
    if specy=='CynomolgusMonkey':
        genes=GeneInfo16.objects.filter(id__iexact=id) 
    if specy=='Sheep':
        genes=GeneInfo17.objects.filter(id__iexact=id)
    if specy=='Rabbit':
        genes=GeneInfo18.objects.filter(id__iexact=id)
    if specy=='Rice':
        genes=GeneInfo19.objects.filter(id__iexact=id)
        #################################################
    if specy=='Rhesusmonkeyhamster':
        genes=GeneInfo20.objects.filter(id__iexact=id)
    if specy=='Bakersyeast':
        genes=GeneInfo21.objects.filter(id__iexact=id)
    if specy=='Fissionyeast':
        genes=GeneInfo22.objects.filter(id__iexact=id) 
    if specy=='Pig':
        genes=GeneInfo23.objects.filter(id__iexact=id)
    if specy=='Breadwheat':
        genes=GeneInfo24.objects.filter(id__iexact=id)
    if specy=='Winegrape':
        genes=GeneInfo25.objects.filter(id__iexact=id)
    if specy=='Westernclawedfrog':
        genes=GeneInfo26.objects.filter(id__iexact=id) 
    if specy=='Maize':
        genes=GeneInfo27.objects.filter(id__iexact=id)



         
    return render_to_response('genedatabase/gene_details.html',{
        'local':local,
        'thisyear':thisyear,
        'genes':genes,
        },context_instance=RequestContext(request))
        
        
        
def loadmore(request,specy):
    if request.method=="POST":
        #异步刷新获取数据
        return HttpResponse("欢迎使用ajax")
        if specy=="Human":
            genes=GeneInfo.objects.all()
            pagesize = int(request.GET.get('ps', '10'))
            paginator = Paginator(genes,pagesize)#使用paginator对象
            page = int(request.GET.get('p', '1'))#取当前页的号码
            genes = paginator.page(page).object_list
            return HttpResponse("欢迎使用ajax")
    return render_to_response('genedatabase/load_more.html',{
        'local':local,
        'thisyear':thisyear,
        'stype':specy
        },context_instance=RequestContext(request))
        
        
        