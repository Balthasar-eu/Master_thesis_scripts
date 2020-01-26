#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 20:42:10 2019

@author: balthasar
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import from_levels_and_colors
from matplotlib import animation
from matplotlib import colors

virulence = {
"ERR1100909" : 50.9,
"ERR1100911" : 71.3,
"ERR1100912" : 37.1,
"ERR1100913" : 61.5,
"ERR1100923" : 41.9,
"ERR1100924" : 58.1,
"ERR1100925" : 40.7,
"ERR1100926" : 31,
"ERR1100928" : 25,
"ERR1100930" : 50,
"ERR1100935" : 52.2,
"ERR1100936" : 65.5,
"ERR1100937" : 46.2,
"ERR1100938" : 65.5,
"ERR1100939" : 71.3,
"ERR1100940" : 45,
"ERR1100941" : 45.2,
"ERR1100942" : 61.5,
"ERR1100944" : 52.2,
"ERR1100946" : 29.3,
"ERR1100947" : 7,
"ERR1100949" : 68.2,
"ERR1100950" : 16,
"ERR1100952" : 18.2,
"ERR1100954" : 61.5,
"ERR1100955" : 40,
"ERR1100972" : 16,
"ERR1100973" : 16,
"ERR1100975" : 16,
"ERR1100976" : 16,
"10CEB535LM" : 7,
"10CEB540LM" : 7,
"10CEB550LM" : 50.9,
"10CEB552LM" : 50.9,
"10CEB553LM" : 45,
"10CEB554LM" : 45,
"10CEB559LM" : 33.3,
"10CEB560LM" : 33.3,
}

clinical_freq = np.array([50.9,
71.3,
37.1,
61.5,
41.9,
58.1,
40.7,
31,
25,
50,
52.2,
65.5,
46.2,
65.5,
71.3,
45,
45.2,
61.5,
52.2,
29.3,
7,
68.2,
16,
18.2,
61.5,
40,
16,
16,
16,
16,
7,
7,
50.9,
50.9,
45,
45,
33.3,
33.3])


arla_class = np.array([2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 3,
       4, 2, 3, 4, 2, 2, 2, 2, 3, 3, 2, 2, 1, 3, 3, 3, 2, 2, 2, 3, 2, 3,
       3, 3, 3, 2, 2, 3])

samplenames = ["ERR1100909",
"ERR1100911",
"ERR1100912",
"ERR1100913",
"ERR1100923",
"ERR1100924",
"ERR1100925",
"ERR1100926",
"ERR1100928",
"ERR1100930",
"ERR1100935",
"ERR1100936",
"ERR1100937",
"ERR1100938",
"ERR1100939",
"ERR1100940",
"ERR1100941",
"ERR1100942",
"ERR1100944",
"ERR1100946",
"ERR1100947",
"ERR1100949",
"ERR1100950",
"ERR1100952",
"ERR1100954",
"ERR1100955",
"ERR1100972",
"ERR1100973",
"ERR1100975",
"ERR1100976",
"10CEB535LM",
"10CEB540LM",
"10CEB550LM",
"10CEB552LM",
"10CEB553LM",
"10CEB554LM",
"10CEB559LM",
"10CEB560LM",]

plt.close("all")

pcarray = np.array([[1.702317345467987131e-01,1.793074845661568051e-01,1.591785352437785905e-01,1.913726076332488146e-01,1.827264015223549198e-01,1.889587937359714365e-01,1.418613926272397352e-01,1.637599926383253301e-01,1.596238915243560486e-01,1.827932317748768853e-01,1.573088855837634081e-01,1.709139140405695756e-01,1.265923001987530605e-01,1.728929842896226310e-01,1.714328176923065883e-01,1.547290753175961653e-01,1.648990091290427895e-01,1.828896403016933847e-01,1.824982656571304818e-01,1.525107527547746589e-01,1.488971515178054705e-01,1.468449826958658500e-01,1.506384982297682695e-01,1.165221429704761252e-01,1.828417877125454449e-01,1.355212917187440802e-01,1.677653242015234858e-01,1.696871837014019557e-01,1.417770270081684103e-01,1.746533640384739516e-01,1.405651851445729184e-01,1.507921496960334118e-01,1.626034807351347533e-01,1.624080543156807832e-01,1.250814503506584296e-01,1.334166441701416517e-01,1.782110369899398972e-01,1.780639245937856718e-01],
[1.979430632961178271e-01,1.579028625311761358e-01,4.361543603146361692e-02,1.652160418324551128e-01,1.055718946676076492e-01,1.760183718065155833e-01,-7.692594137156266720e-02,-1.153730091535467989e-01,-1.806921005315596240e-01,1.553707508071759114e-01,3.358164728082218414e-02,1.647982822303895867e-01,-1.393705711181285334e-01,1.649373441832392428e-01,1.509165826260730481e-01,-1.865433887070183649e-01,2.085442346733507729e-02,1.596512283699320633e-01,6.946620096242708819e-02,-1.112843779698753366e-01,-2.051695256818253932e-01,-1.554103013907791320e-01,-2.131161985088519673e-01,-1.476121877597346144e-01,1.596928751925159873e-01,-7.107776273229861330e-02,-2.051057000178916345e-01,-1.547342622356585040e-01,-2.156556935469362279e-01,-2.145215317869174354e-01,-2.116540461692386355e-01,-1.499074016639973250e-01,1.856110106461009435e-01,1.853555365828266910e-01,-2.099684057477160248e-01,-2.074806534497337651e-01,1.832970674606111128e-01,1.831972700956197853e-01],
[-2.193328622959289098e-01,5.116624827352262567e-02,1.689477016056521175e-01,6.581651984539044209e-02,9.088064728725091246e-02,-6.850278136635240024e-03,-3.134676426337664606e-02,1.481824953186143545e-02,-2.131203893907428551e-02,3.419046151522640159e-03,1.821142293313008156e-02,7.915990286555060762e-02,-6.672647632406619733e-02,7.791082077117414406e-02,8.805279289559997391e-02,-2.397336253525696875e-03,1.766028675162261974e-01,8.931372274198531858e-02,1.015836753822769029e-01,-1.601074442222704952e-01,-3.096942964765020401e-01,-5.158183739425128284e-02,2.279009502094343775e-01,-3.553615549555953768e-02,8.928314674834401943e-02,2.332638401519680424e-02,1.937567675445876492e-01,2.395238458518026969e-01,2.157587149484504474e-01,1.903732009836238692e-01,-3.768993851231546288e-01,-3.595063665985388512e-01,-3.072717281636825892e-01,-3.064762158964081107e-01,-3.336874150202651101e-02,-8.053614625744758004e-02,-2.481508149801994398e-02,-2.396807012524968666e-02]])

arla_pca = np.array([
[8.957519064885077764e-02,9.517746487229712671e-02,9.432054183546959014e-02,9.858743243974932291e-02,1.079271930570622345e-01,9.962544839807536967e-02,8.327227254161126413e-02,9.173537976788781467e-02,9.408882028354437210e-02,9.895831257130939529e-02,1.091864286943635826e-01,8.797373092953741924e-02,7.514999545791686830e-02,8.913954092689291286e-02,8.647110074181076078e-02,8.578886265648902332e-02,9.637529885373453953e-02,9.103364033185232240e-02,1.081246056618085927e-01,8.999253680830564261e-02,9.548007340178624180e-02,9.415021489678372657e-02,9.402448623346694456e-02,6.509048783993630349e-02,9.101405407739733555e-02,7.345434862756471117e-02,9.664027891899659717e-02,9.240370708390474586e-02,9.243999245732083259e-02,9.646025281666596729e-02,9.695851224288899728e-02,1.031517032576402260e-01,8.856538871116582945e-02,8.822789827134031038e-02,7.387961525514161187e-02,8.064912239971257224e-02,9.828235672995087879e-02,9.798377419797797694e-02,9.960168181914738628e-02,1.310318574348513510e-01,9.964415057403207265e-02,1.310259415747268563e-01,1.292374370039547915e-01,1.306198118318253609e-01,1.291602282930560397e-01,9.994466478001884369e-02,8.428160209544882653e-02,1.378258304552981428e-01,1.392553883297030504e-01,9.714792234722058328e-02,1.049762084271670531e-01,1.039216069067789849e-01,1.032073071307653273e-01,9.701906236637242165e-02,1.038287704242474158e-01,1.001542962859332531e-01,1.050625316203281229e-01,9.887006236632819856e-02,8.412323383764536144e-02,1.390342922857845387e-01,9.866716672970360369e-02,3.580778290177635131e-02,1.392536898145836366e-01,1.021904301679020660e-01,9.598201150601989862e-02,1.013059879278933040e-01,1.011436264772079957e-01,9.598263524429614946e-02,1.393517538152836688e-01,1.392736166862242331e-01,9.703319869898271999e-02,9.700359979883119155e-02,1.065075592813587851e-01,1.393341645436293741e-01,1.392593969482082139e-01,1.392493349696610205e-01,1.000446265415482344e-01,9.315881710040020691e-02,9.525779728003989921e-02,1.355803897411097525e-01,1.377269393217154525e-01,1.357871102995873991e-01,1.357825846069663778e-01,1.358122575042801650e-01,1.357664035171446959e-01,1.377358229850838001e-01,1.107894047140914190e-01,1.392486380338435192e-01],
[-1.264806537873973125e-01,-8.206925594963652981e-02,-2.826080872761396745e-02,-8.776648251557717106e-02,-6.231703261840308011e-02,-9.787068659258560399e-02,8.828040047843381533e-02,1.108517788646030400e-01,1.435575049510319223e-01,-8.852787758274932584e-02,-7.421383414073992602e-02,-1.006835228558975309e-01,1.104378374588794315e-01,-1.006439773556218747e-01,-7.854070029517051799e-02,1.482906873267054171e-01,-2.536304211219870750e-02,-8.594228736784928491e-02,-5.459362816745982894e-02,8.865492542410362276e-02,1.126480285528519143e-01,1.152607570721539132e-01,1.340552876932504800e-01,1.187847014676304824e-01,-8.599010623362710048e-02,8.218093667960355009e-02,1.401519676096394662e-01,1.297914407993404906e-01,1.210780971728001310e-01,1.555349825652648510e-01,1.170661341412086442e-01,8.224969957803410225e-02,-1.233372896361680515e-01,-1.232734415142285184e-01,1.620338472531331342e-01,1.564564069163902893e-01,-1.070243074523883758e-01,-1.071004625838358898e-01,1.189412234681617303e-01,-8.125329285654792066e-02,1.189649853989028228e-01,-8.122281538931259981e-02,-4.905235097077387252e-02,-4.691397564423743088e-02,-4.899466219239463627e-02,1.180751188550243497e-01,8.941882408143778105e-02,-9.819493193471352799e-02,-9.623549988907650898e-02,1.358562363290541364e-01,8.741669647156746770e-02,-9.321231819606680402e-02,-2.269951090453104989e-02,1.468795255727415872e-01,1.257268933467942384e-01,1.182664405949096359e-01,1.400548914641486509e-01,-1.057522337928099437e-01,8.957497410340214961e-02,-9.661746528847076076e-02,-1.057427482826956699e-01,1.698512861520167644e-02,-9.624457485407544866e-02,-1.043387234626107635e-01,1.393611684419897856e-01,1.286191048152315886e-01,1.287733687338246669e-01,1.393611515838752424e-01,-9.651762996110903325e-02,-9.640133887661980761e-02,1.239896099470401691e-01,1.239227062203308721e-01,1.306136435303284771e-01,-9.655108257585978682e-02,-9.626063914559250534e-02,-9.623773009363388586e-02,1.041976721122148120e-01,1.368223453519144261e-01,1.351831558881613216e-01,-9.603826439589413511e-02,-9.883514877597172898e-02,-9.659078036815837986e-02,-9.641336283433450671e-02,-9.659837793569157893e-02,-9.648647641214605497e-02,-9.883112982155861648e-02,6.742925027202019417e-02,-9.623838425323062717e-02],
[1.590435348246941050e-01,1.755487919027874932e-01,7.880069546751862497e-02,1.850093721969744143e-01,1.121659613400545485e-01,1.839480538667164766e-01,7.640824099998913543e-02,8.767372146112728526e-02,3.909403390756251845e-02,1.547188611964382676e-01,-6.214483745893822908e-02,1.931116488146789600e-01,1.527646063935757595e-02,1.944902279834091197e-01,1.783982422976212734e-01,4.136747583892283475e-02,6.208339206556356960e-02,1.828679288398155356e-01,7.952062592027156485e-02,3.405865126561656575e-02,-6.954432879347034746e-02,-3.995258627024854577e-03,-3.352957865421214712e-02,3.465172581449654787e-02,1.827905698565285564e-01,8.924510623266734832e-02,2.018957178356889023e-02,8.161180961314710991e-02,-6.637501142333861182e-02,4.550610581077654121e-02,-9.365479334261729571e-02,-6.286644717893441481e-02,1.359738324706946944e-01,1.360808643557386288e-01,9.166495628391051198e-03,-5.508717380307404272e-03,1.576262762169381992e-01,1.577312837232106302e-01,3.035465364950283573e-02,-1.790697361522084474e-01,3.040029870412800511e-02,-1.790699713067399890e-01,-1.620343282412589625e-01,-1.840006055495959358e-01,-1.620343444009101075e-01,-9.474639441142027399e-03,5.977609752587931258e-02,-1.202057828868767081e-01,-1.160661686878403931e-01,4.927068134467142008e-02,-4.248302871331453845e-03,1.756423363658207448e-01,2.307006419732041125e-02,-1.479833770726794208e-02,-1.721436111194266788e-02,2.219893803510691183e-02,3.783503217601239671e-03,1.808130493594477528e-01,5.982257977976172120e-02,-1.158013515629484624e-01,1.804977088477403291e-01,-5.499820491738294931e-03,-1.160922615888282244e-01,1.299106146615431057e-01,3.118943554709899257e-02,6.090995612927896895e-03,6.145290680212790041e-03,3.119026381612189855e-02,-1.159423697212575571e-01,-1.165538471418502003e-01,4.642643433351489185e-02,4.629443379004812703e-02,1.016682517269357731e-02,-1.156994865097135811e-01,-1.160866427423581093e-01,-1.160874298011638772e-01,1.103951907687477155e-02,-3.008520653778270701e-02,2.564940207936970695e-02,-9.865516560236715915e-02,-1.205011056187877128e-01,-9.843183102006831942e-02,-9.859310959023968868e-02,-9.837496718190985590e-02,-9.837812778581732365e-02,-1.204409125322809943e-01,3.072192664080875685e-02,-1.160848046824676377e-01]])

newv = np.array([(1 if virulence[v]<25 else (2 if virulence[v] < 50 else (3 if virulence[v] < 60 else 4)) ) for v in virulence])

from matplotlib.cm import plasma
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.ioff()
animate = False

def rotate(i):
    ax.view_init(elev=10., azim=i)
    return fig,

ffmpegwriter = animation.FFMpegWriter(fps=30, codec="theora", extra_args=["-qscale:v", "10"])#, codec="theora",bitrate=5000 ) #, codec=None, bitrate=None, extra_args=None, metadata=None)

###############################################################################
####################### Histogram freuqency distribution ######################

fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)
ax.hist(clinical_freq, np.arange(0,105,5))
#fig.colorbar(plasma)

ax.set_xlabel("clinical frequency")
ax.set_ylabel("count")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels)
fig.savefig("plots/histogram_freq.pdf")

###############################################################################
######################### Histogram class distribution ########################

fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)
#ax.plot (clinical_freq)
#ax.hist([clinical_freq[newv==1],clinical_freq[newv==2],clinical_freq[newv==3],clinical_freq[newv==4]],np.arange(0,105,5), label="class")

rangedict = {
        1: "  0 - 25",
        2: " 25 - 50",
        3: " 50 - 60",
        4: " 60 -100"
        }

[ax.hist(clinical_freq[newv==i] ,np.arange(0,105,5), label="class "+str(i)+rangedict[i]) for i in range(1,5)]

ax.set_xlabel("clinical frequency")
ax.set_ylabel("count")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels)
fig.savefig("plots/histogram.pdf")

###############################################################################
######################### Histogram class distribution ########################

fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)

classcount = [sum(newv==i) for i in range(1,5)]

for i, c in enumerate(classcount):
    temp = [0,0,0,0]
    temp[i] = c
    ax.bar([1,2,3,4], temp, 0.5) 

ax.set_xlabel("class")
ax.set_ylabel("count")

ax.set_xticks([1,2,3,4])
ax.set_xticklabels(["1","2","3","4"])

for x,y in enumerate(classcount,1):
    ax.text(x,y/2, y,
            ha='center', va='center', color = 'white',
            path_effects = [path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

fig.savefig("plots/count_histogram.pdf")


###############################################################################
################################# PCA our data ################################
################################### classes ###################################

fig = plt.figure(figsize=(5, 5), frameon=False)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
[ax.scatter(pcarray[0, newv==c],pcarray[1, newv==c],pcarray[2, newv==c], label="class "+str(c)) for c in [1,2,3,4]] 
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels)

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

for x,y,z,n in zip(pcarray[0],pcarray[1],pcarray[2],newv):
    ax.text(x,y,z,n)

fname = "pca_class"

fig.savefig('plots/' + fname + '.pdf')

if animate:
    with ffmpegwriter.saving(fig,'plots/' + fname + '.avi', 300):
        for angle in range(0, 360):
            ax.view_init(30, angle)
            ffmpegwriter.grab_frame()
    
###############################################################################
################################# PCA our data ################################
############################### clinical freqency #############################

fig = plt.figure(figsize=(5, 5), frameon=False)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
#ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(pcarray[0],pcarray[1],pcarray[2], c=clinical_freq,cmap=plasma, label="Clinical frequency")
#fig.colorbar(p, cax=ax)

axins1 = inset_axes(ax,
                    width="50%",  # width = 50% of parent_bbox width
                    height="5%",  # height : 5%
                    loc='upper right')

fig.colorbar(p, cax=axins1, orientation="horizontal")
axins1.xaxis.set_ticks_position("bottom")


ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

for x,y,z,n in zip(pcarray[0],pcarray[1],pcarray[2],clinical_freq):
    ax.text(x,y,z,n)
    
    
fname = "pca_freq"

fig.savefig('plots/' + fname + '.pdf')

if animate:

    with ffmpegwriter.saving(fig,'plots/' + fname + '.avi', 300):
        for angle in range(0, 360):
            ax.view_init(30, angle)
            ffmpegwriter.grab_frame()

###############################################################################
########################### PCA our data plus Arla ############################

fig = plt.figure(figsize=(5, 5), frameon=False)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
cclass = np.concatenate((newv,arla_class))
#ax.scatter(arla_pca[0],arla_pca[1],arla_pca[2], c = np.concatenate((newv,arla_class)))
[ax.scatter(arla_pca[0, cclass==c],arla_pca[1, cclass==c],arla_pca[2, cclass==c], label="class "+str(c), color=color) for color,c in zip(["blue","green","orange","red"],[1,2,3,4])]
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels)

fig.patch.set_visible(False)
ax.patch.set_visible(False)

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

for x,y,z,n in zip(arla_pca[0],arla_pca[1],arla_pca[2],list(newv)+["A"]*50):
    ax.text(x,y,z,n)

fname = "pca_arla"

fig.savefig('plots/' + fname + '.pdf')
fig.savefig('plots/' + fname + '.png')

if animate:
    with ffmpegwriter.saving(fig,'plots/' + fname + '.avi', 300):
        for angle in range(0, 360):
            ax.view_init(30, angle)
            ffmpegwriter.grab_frame()

plt.close("all")

###############################################################################
########################### PCA Arla only          ############################

plt.ion()

newno = np.array([99 for v in virulence])

fig = plt.figure(figsize=(5, 5), frameon=False)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
cclass = np.concatenate((newno,arla_class))
#ax.scatter(arla_pca[0],arla_pca[1],arla_pca[2], c = np.concatenate((newv,arla_class)))
[ax.scatter(arla_pca[0, cclass==c],arla_pca[1, cclass==c],arla_pca[2, cclass==c], label="class "+str(c), color=color) for color,c in zip(["blue","green","orange","red"],[1,2,3,4])]
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels)

fig.patch.set_visible(False)
ax.patch.set_visible(False)

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

for x,y,z,n in zip(arla_pca[0],arla_pca[1],arla_pca[2],[""]*len(newv)+list(range(1,51))):
    ax.text(x,y,z,n)

fname = "pca_arla_only"

fig.savefig('plots/' + fname + '.pdf')
fig.savefig('plots/' + fname + '.png')

if animate:
    with ffmpegwriter.saving(fig,'plots/' + fname + '.avi', 300):
        for angle in range(0, 360):
            ax.view_init(30, angle)
            ffmpegwriter.grab_frame()

#plt.close("all")

###############################################################################
########################### Arla histogram ############################


Arla_results = np.array([2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 3, 4, 2, 3, 4, 2, 2, 2, 2, 3, 3, 2, 2, 1, 3, 3,
 3, 2, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 3])



fig = plt.figure(frameon=False)
ax = fig.add_subplot(111)

classcount = [sum(Arla_results==i) for i in range(1,5)]
color = ["blue", "green", "orange", "red"]

for i, c in enumerate(classcount):
    temp = [0,0,0,0]
    temp[i] = c
    ax.bar([1,2,3,4], temp, 0.5, color=color[i]) 

ax.set_xlabel("class")
ax.set_ylabel("count")

ax.set_xticks([1,2,3,4])
ax.set_xticklabels(["1","2","3","4"])

for x,y in enumerate(classcount,1):
    ax.text(x,y/2, y,
            ha='center', va='center', color = 'white',
            path_effects = [path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

fig.savefig("plots/histogram_arla.svg")
fig.savefig("plots/histogram_arla.pdf")


#fig = plt.figure(frameon=False)
#ax = fig.add_subplot(111)
#ax.hist(Arla_results)
##fig.colorbar(plasma)
#
#ax.set_xlabel("class")
#ax.set_ylabel("count")

plt.ioff()
