import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import matplotlib.ticker as mtick


#################################################################

dates = pd.date_range(start="3/4/2020",end="5/10/2020",freq='D')

 
# Values of each group
bars1 = [3,0,2,2,2,3,1,3,3,5,6,6,10,6,15,11,15,28,32,18,33,35,38,37,62,31,39,
            25,49,30,41,39,10,49,65,76,174,106,94,33,40,52,59,92,61,72,57,102,
            62,106,91,48,52,74,62,70,42,82,77,52,35,30,42,36,26,35,49,20]
bars2 = [0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,1,2,0,3,3,5,2,1,5,4,7,6,8,11,8,13,
              14,4,25,13,9,35,13,6,14,15,15,14,18,10,9,10,11,8,14,9,11,5,9,6,8,
              6,6,2,4,2,0,4,3,2,0,1,1]
 
# Heights of bars1 + bars2
bars = np.add(bars1, bars2).tolist()
 
# The position of the bars on the x-axis
r = np.arange(68)
 
# Names of group and bar width
names = dates.date
barWidth = 0.9

plt.ylim(top=230)
 
p1 = plt.bar(r, bars2, color='pink', width=barWidth,zorder=2)
p2 = plt.bar(r, bars1, bottom=bars2, color='gray', width=barWidth,zorder=2)

# Custom X axis
plt.xticks(r[r%2==0], names[r%2==0], rotation=90, fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel("number of cases", fontsize=22)

# plt.locator_params(axis='x', nbins=34)
plt.grid(axis='y',linestyle=':',zorder=5)

fig = plt.gcf()
fig.set_size_inches(20, 10)

ax = plt.gca()
rects = ax.patches

for rect, label, hp in zip(rects, bars2, bars1):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height + hp + 1, label,
            ha='center', va='bottom', color='red', fontsize=13)
    
for rect, label in zip(rects, bars1):
    height = rect.get_height()
    ax.text(rect.get_x() + rect.get_width() / 2, height + label + 7, label,
            ha='center', va='bottom', fontsize=13)

plt.legend((p1[0], p2[0]), ('deaths', 'confirmed COVID-19 cases'), fontsize=22)
 
# Show graphic
plt.savefig('epicurve.pdf', format="pdf", bbox_inches='tight')
plt.close()



#################################################################



bp = [0.238753230629652,0.298441538287065,0.417818153601891,0.59688307657413,
      0.59688307657413,0.656571384231543,0.716259691888956,0.895324614861195,
      1.07438953783343,1.13407784549085,1.43251938377791,1.67127261440756,
      1.91002584503722,2.26815569098169,2.80535045989841,3.34254522881513,
      4.1184932283615,4.29755815133374,5.43163599682458,5.9688307657413,
      7.10290861123215,8.47573968735265,10.0276356864454,10.8632719936492,
      11.2214018395936,11.9376615314826,12.2361030697697,12.9523627616586,
      14.0267522994921,14.6833236837236,14.1461289148069,15.45927168327,
      16.8321027593905,18.6824402967703,26.740361830521,29.7247772133917,
      29.3069590597898,29.9038421363639,31.0379199818548,32.3510627503178,
      33.7835821340958,37.3648805935405,38.3198935160591,39.5136596692074,
      40.9461790529853,43.9902827435134,43.5724645899115,44.2887242818004,
      37.9020753624572,35.1564132102162,35.4548547485033,34.6789067489569,
      34.0820236723828,34.3804652106699,33.5448289034661,31.3960498277992,
      31.0976082895122,30.1425953669936,29.0085175215027,25.4272190620579,
      24.8303359854838,22.6218686021595,20.234336295863,18.9808818350573,
      18.3839987584832,18.1452455278536]

cs = [0.0121330157071595,0.024266031414319,0.024266031414319,0.024266031414319,
      0.0363990471214784,0.0606650785357974,0.0849311099501163,
      0.121330157071595,0.169862219900233,0.242660314143189,0.315458408386146,
      0.363990471214784,0.46105459687206,0.533852691115017,0.606650785357974,
      0.812912052379685,1.05557236652287,1.23756760213027,1.45596188485914,
      1.75928727753812,1.96554854455983,2.12327774875291,2.52366726708917,
      2.72992853411088,3.06965297391135,3.22738217810442,3.7127028063908,
      3.88256502629103,4.11309232472706,4.28295454462729,4.13735835614138,
      4.51348184306332,4.71974311008504,4.92600437710675,5.35065992685733,
      5.67825135095063,6.17570499494417,6.16357197923701,6.05437483787258,
      6.19997102635849,6.06650785357974,6.21210404206565,6.22423705777281,
      6.32130118343009,6.67315863893771,6.5275624504518,6.51542943474464,
      6.79448879600931,6.77022276459499,6.60036054469475,6.0179757907511,
      6.61249356040191,6.89155292166658,6.96435101590954,6.83088784313078,
      7.00075006303102,7.15847926722409,7.04928212585965,6.9158189530809,
      6.63675959181623,6.466897371916,5.93304468080098,5.54478817817188,
      5.5083891310504,5.54478817817188,4.84107326715663]

plt.yscale("log")
plt.ylim(top=100)
plt.ylim(bottom=0.01)
plt.plot(range(0,len(bp),1), bp,linewidth=3.0)
plt.plot(range(0,len(bp),1), cs,linewidth=3.0)
plt.legend(('Budapest', 'Other regions'), fontsize=22)
plt.grid(axis='y',linestyle=':',zorder=5)
plt.ylabel("14-day cumulative incidence (per\n100,000 population), log scale", fontsize=22)

dates = pd.date_range(start="3/6/2020",end="5/10/2020",freq='D')

r = np.arange(len(bp))
names = dates.date

plt.xticks(r[r%2==0], names[r%2==0], rotation=90, fontsize=20)
plt.yticks([0.01, 0.1, 1, 10, 100], [0.01, 0.1, 1, 10, 100], fontsize=20)

fig = plt.gcf()
fig.set_size_inches(20, 10)
plt.savefig('14day_cum_inc.pdf', format="pdf", bbox_inches='tight')
plt.close()


#######################################################################


hely = [0, 0, 0, 1.12, 1.74, 4.87, 20.85, 28.57]
hol = range(len(hely))
also = [0, 0, 0, 1, 1, 1, 3, 3]
felso = [0, 0, 0, 2, 2, 2, 3, 3]
meret = np.array([0, 0, 0, 3, 8, 35, 177, 198])*10
xcimke = ["<10 years", "10-19", "20-29", "30-39", "40-49", "50-64", "65-79", "80+"]

fig = plt.gcf()
fig.set_size_inches(15, 10)

plt.grid(axis='y',linestyle=':',zorder=5)
plt.scatter(hol, hely, s=meret,zorder=3, color='orange')
plt.errorbar(hol, hely, yerr=[also,felso], color='orange', fmt='.', markersize=0, ecolor='black',capsize=4, elinewidth=1)
plt.xticks(hol, xcimke, fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel("crude CFR", fontsize=20)

ax = plt.gca()
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))

plt.savefig('crudeCFR_agegroup.pdf', format="pdf", bbox_inches='tight')

