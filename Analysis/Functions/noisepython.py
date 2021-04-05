import csv
import seaborn as sns
import pandas
import matplotlib.pyplot as plt

sns.set(style="whitegrid")

noiseintegrated = []
noisemax = []
noisetime = []
noisecat = []
rownum = 1
with open('C:\\Users\\Parthiv_Lab\\Documents\\Data\\180312_Cos_Noise\\noise_nonoisestats.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        if rownum >1:
            if float(row[1]) >= .05 and float(row[0])<100 and float(row[2]) < 120:
                noiseintegrated.append(float(row[0]))
                noisemax.append(float(row[1]))
                noisecat.append(row[3])
                noisetime.append(float(row[2]))
#                if float(row[0]) >= 10:
 #                   noisetime.append(float(row[2]))
  #              else:
   #                 noisetime.append(float(0))

        rownum += 1

d = {'Integrated Nuclear Fluorescence':noiseintegrated,'Peak Height':noisemax,'Peak Time':noisetime, "Stimulation":noisecat}
df = pandas.DataFrame(data = d)

g = sns.catplot(x = 'Stimulation',y = 'Peak Time',kind="box", data=df)
sns.stripplot(x = 'Stimulation',y = 'Peak Time',color="k", size=3, data=df, ax=g.ax, jitter=.3);

plt.savefig('C:\\Users\\Parthiv_Lab\\Documents\\Data\\180312_Cos_Noise\\PeakTime_1.svg', format='svg')
plt.show()
