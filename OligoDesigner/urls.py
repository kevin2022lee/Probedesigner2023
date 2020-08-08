from django.conf.urls import *
from views import *
from geneviews import *
from orderviews import *
from SightViews import *
from blastviews import *
from stxviews import *
from movedesignviews import *
from django.contrib import admin
import settings
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^site_medias/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT }),
    url(r'^help/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.HELP_ROOT }),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^$', index),
    url(r'^eindex', eindex),
    url(r'^probedesign', probedesign),
    url(r'^subseq', subseq),
    url(r'^intro', intro),
#    url(r'^entrez', entrez),
    url(r'^fromfile', fromfile),
    url(r'^english', english),
    url(r'^downloadentrez', downloadentrez),
    url(r'^xmltoentrez', entreztoxml),
    url(r'^seqidtoentrez', entrezseqidtoxml),
    url(r'^convertdata', convertdata),
    url(r'^test', test),
    url(r'^x4merCalc', x4merCalc),
    url(r'^startdesign', startdesign),
    url(r'^bingodesign', bingodesign),
    url(r'^NonNshFilter', NonNshFilter),
    url(r'^Probelist', Probelist),
    url(r'^PostCalcXmer', PostCalcXmer),
    url(r'^Universalvalue', Universalvalue),
    url(r'^ProbeSetsXmer',ProbeSetsXmer),
    url(r'^Probegroupsvalues',Probegroupsvalues),
    url(r'^GenerateProbesets',GenerateProbesets),
    url(r'^GenerateProbeset',GenerateProbeset),
    url(r'^Probesetsgenerate',Probesetsgenerate),
    url(r'^bigdata',BigData),
    url(r'^tfdeal',tfdeal),
    url(r'^datatfshow',tfdealdata),
    url(r'^checkAccessionLen/(.+)/$',checkAccessionLen),
    url(r'^updatelog',updatelog),
    url(r'^sightintro',sightintro),
    url(r'^zzprobe',zzprobe),
    url(r'^seqresolve',seqresolve),
    url(r'^entrezremote',entrezremote),
    url(r'startzzprobedesign', startzzprobedesign),
    url(r'NNFzzprobe',NNFzzprobe),
    url(r'PCXzzprobe',PCXzzprobe),
    url(r'zzProbeSetsXmer',zzProbeSetsXmer),
    url(r'GeneratezzProbesets',GeneratezzProbesets),
    url(r'designbranch',designbranch),
    url(r'branchscorecalc',branchscorecalc),
    url(r'showbranchscore',showbranchscore),
    url(r'genesearch/$',genesearch),
    url(r'gdbsearch/',gdbsearch),
    url(r'sightrnasearch/',sightrnasearch),
    url(r'genedatabase/(\w+)/(\d+)/',geneshow),
    url(r'sightrnadatabase/(\w+)/(\d+)/',sgeneshow),
    url(r'loadmore/(\w+)/',loadmore),
    url(r'tools/msa/',MSA),
    url(r'tools/msaresult/',multiseqalign),
    url(r'order/',order),
    url(r'16xpd/',stxpd),
    url(r'fileparse/',fileparse),
    url(r'entrezparse/',entrezparse),
    url(r'16xpd_start/',stxpd_design),
    url(r'stxprobefilter/',stxprobefilter),
    url(r'stxprobeXmers/',stxprobeXmers),
    url(r'sightRNAsearch/$',sightRNAsearch),
    url(r'gdbsearch/',gdbsearch),
    url(r'move_design/',move_design)
)
